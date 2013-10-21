/*
 * Wenger CD, Coon JJ. A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, Journal of Proteome Research, 2013; 12(3): 1377-86
 * http://www.chem.wisc.edu/~coon/software.php#morpheus
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 * Altered by Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 */
using System;
using System.Collections.Generic;
using System.Xml;
using Proteomics.Utilities;

namespace Trinity
{
    /// <summary>
    /// Loads and stores the sample information sored in a Raw file
    /// </summary>
    public partial class Spectra : List<ProductSpectrum>
    {
        public double MaxMZ;
        public double MinMZ;
        public double MaxProductMass;
        public double MinProductMass;
        public double MaxRt;
        public double MinRt;

        private const bool HARMONIC_CHARGE_DETECTION = false;
        public int NbScans;

        public Tracks tracks;
        public List<MS1Spectrum> MS1s;

        public Spectra(int nbScans = -1) : base() 
        {
            this.NbScans = nbScans;
            MaxMZ = double.MinValue;
            MinMZ = double.MaxValue;
            MaxProductMass = double.MinValue;
            MinProductMass = double.MaxValue;
            MaxRt = double.MinValue;
            MinRt = double.MaxValue;
            tracks = new Tracks();
            MS1s = new List<MS1Spectrum>();
        }

        public void ExportTracks(string filename)
        {
            tracks.Export(filename);
        }

        public void ExportMSMS(string filename)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine(ProductSpectrum.TITLE);

            foreach (ProductSpectrum spectrum in this)
                writer.AddLine(spectrum.ToString());
            writer.writeToFile();
        }

        public static Spectra Import(string filenameMSMS, string filenameTracks)
        {
            Spectra spectra = new Spectra();
            vsCSV csv = new vsCSV(filenameMSMS);
            if (csv.LINES_LIST.Count == 0 || csv.LINES_LIST[0].CompareTo(ProductSpectrum.TITLE) != 0)
                return null;
            for (int i = 1; i < csv.LINES_LIST.Count; i++)
            {
                string[] splits = csv.LINES_LIST[i].Split(vsCSV._Generic_Separator);
                double mz = double.Parse(splits[3]);
                int charge = int.Parse(splits[5]);
                int nbPeaks = int.Parse(splits[9]);
                GraphML_List<MsMsPeak> peaks = new GraphML_List<MsMsPeak>(nbPeaks);
                i++;
                for(int j = 0; j < nbPeaks; i++,j++)
                {
                    try
                    {
                        string[] splitPeaks = csv.LINES_LIST[i].Split('\t');
                        if (splitPeaks.Length > 2)
                            peaks.Add(new MsMsPeak(double.Parse(splitPeaks[0]), double.Parse(splitPeaks[1]), int.Parse(splitPeaks[2])));
                        else
                            peaks.Add(new MsMsPeak(double.Parse(splitPeaks[0]), double.Parse(splitPeaks[1]), 0));
                    }
                    catch (Exception)
                    {
                        Console.WriteLine("Error parsing line : " + csv.LINES_LIST[i]);
                    }
                }
                spectra.AddMSMS(new ProductSpectrum(int.Parse(splits[0]), double.Parse(splits[1]), splits[2], mz, double.Parse(splits[4]), charge,Proteomics.Utilities.Numerics.MassFromMZ(mz, charge), peaks, double.Parse(splits[8])));                
            }
            if(!string.IsNullOrEmpty(filenameTracks))
                spectra.tracks = Tracks.Import(filenameTracks);
            return spectra;
        }

        public void AddMSMS(ProductSpectrum spectrum)
        {
            if (spectrum.PrecursorMZ > MaxMZ)
                MaxMZ = spectrum.PrecursorMZ;
            if (spectrum.PrecursorMZ < MinMZ)
                MinMZ = spectrum.PrecursorMZ;
            if (spectrum.RetentionTimeInMin > MaxRt)
                MaxRt = spectrum.RetentionTimeInMin;
            if (spectrum.RetentionTimeInMin < MinRt)
                MinRt = spectrum.RetentionTimeInMin;

            foreach (MsMsPeak peak in spectrum.Peaks)
            {
                if (peak.MZ > MaxProductMass)
                    MaxProductMass = peak.MZ;
                if (peak.MZ < MinProductMass)
                    MinProductMass = peak.MZ;
            }
            Add(spectrum);
        }

        public static Spectra Load(pwiz.CLI.msdata.MSDataFile msFile, DBOptions options, string filePath, bool loadMS = true, bool filterMS2 = true)
        {
            //Find file name in msFile;
            string mzMlFilepath = filePath;
            int num_spectra = msFile.run.spectrumList.size();
            Spectra spectra = new Spectra(num_spectra);
            //List<Trail> trails = new List<Trail>();       
            MS1Spectrum previousMS1 = null;
            try
            {

                //TODO DONT forget to remove the limiter
                //int maxNbMSMS = 10;
                for (int i = 0; i < num_spectra/* && i < 200*/; i++)//TODO Fix that later!
                {
                    //Spectrum
                    pwiz.CLI.msdata.Spectrum spec = msFile.run.spectrumList.spectrum(i, true);

                    if (spec.precursors.Count > 0 || spec.cvParam(pwiz.CLI.cv.CVID.MS_ms_level).value > 1)//is an MSMS
                    {
                        double retention_time = spec.scanList.scans[0].cvParam(pwiz.CLI.cv.CVID.MS_scan_start_time).timeInSeconds() / 60.0;

                        //List precursors and their intensities
                        double precursor_mz = 0;//Is there a value for the time a scan took to complete?
                        int charge = 2;
                        double precursor_intensity = 0;
                        string fragmentation_method = "unknown";
                        double isolationWindow = 1.0;
                        foreach (pwiz.CLI.msdata.Precursor precursor in spec.precursors)
                        {
                            fragmentation_method = precursor.activation.cvParams[0].name;
                            if (precursor.isolationWindow.cvParams.Count > 2 && (double)precursor.isolationWindow.cvParams[1].value == (double)precursor.isolationWindow.cvParams[2].value)
                                isolationWindow = precursor.isolationWindow.cvParams[1].value;
                            else if (precursor.isolationWindow.cvParams.Count > 2)
                                Console.WriteLine("Weird Isolation Window");

                            foreach (pwiz.CLI.msdata.SelectedIon ion in precursor.selectedIons)
                            {
                                //Cycle through MS to get real precursor intensities
                                precursor_mz = ion.cvParams[0].value;
                                if (ion.cvParams.Count > 1)
                                    charge = (int)ion.cvParams[1].value;
                                //else
                                //    Console.WriteLine("No charge computed for precursor ");
                                if (ion.cvParams.Count > 2)
                                    precursor_intensity = ion.cvParams[2].value;
                            }
                        }


                        int scan_index = i;
                        int scan_number = scan_index + 1;

                        pwiz.CLI.msdata.BinaryDataArray mz = spec.getMZArray();
                        pwiz.CLI.msdata.BinaryDataArray intensity = spec.getIntensityArray();

                        int num_peaks = mz.data.Count;
                        if (num_peaks != intensity.data.Count)
                        {
                            Console.WriteLine("PreoteWizard reports peaks arrays (mz/intensity) of different sizes : (" + num_peaks + "/" + intensity.data.Count + ")");
                            if (intensity.data.Count < num_peaks)
                                num_peaks = intensity.data.Count;
                        }
                        GraphML_List<MsMsPeak> peaks = new GraphML_List<MsMsPeak>(num_peaks);
                        for (int k = 0; k < num_peaks; k++)
                        {
                            if (intensity.data[k] > 0)
                            {
                                MsMsPeak peak = new MsMsPeak(mz.data[k], intensity.data[k], 0);
                                peaks.Add(peak);
                            }
                        }
                        mz.Dispose(); mz = null;
                        intensity.Dispose(); intensity = null;

                        peaks.Sort(MsMsPeak.AscendingMzComparison);

                        if (filterMS2)
                        {
                            //peaks = AssignChargeStates(peaks, options.maximumAssumedPrecursorChargeState, options.precursorMassTolerance);
                            //peaks = Deisotopebkp(peaks, options.maximumAssumedPrecursorChargeState, options.precursorMassTolerance);
                            peaks = AssignChargeStatesAndDeisotope(peaks, options.MaximumPrecursorChargeState, new MassTolerance(options.productMassTolerance.Value * 0.5, options.productMassTolerance.Units));
                            peaks = FilterPeaks(peaks, options.MaximumNumberOfFragmentsPerSpectrum);

                            //TODO Add Contaminant removal 
                            //peaks = ContaminantMasses.RemoveContaminantsFromMzSortedList(peaks, options.productMassTolerance);

                            //Can sometime be sorted by intensity after this call
                            //peaks = FilterPeaksV2(peaks);
                            peaks.Sort(MsMsPeak.AscendingMzComparison);
                        }

                        /*//TODO Validate that in most cases, next steps can calculate missing charge
                        if (charge == 0)
                        {
                            for (int c = options.minimumAssumedPrecursorChargeState; c <= options.maximumAssumedPrecursorChargeState; c++)
                            {
                                if (options.assignChargeStates)
                                {
                                    peaks = AssignChargeStates(peaks, c, options.productMassTolerance);
                                    if (options.deisotope)
                                    {
                                        peaks = Deisotope(peaks, c, options.productMassTolerance);
                                    }
                                }

                                double precursor_mass = Utilities.MassFromMZ(precursor_mz, c);

                                ProductSpectrum spectrum = new ProductSpectrum(mzMlFilepath, scan_number, retention_time, fragmentation_method, precursor_mz, precursor_intensity, c, precursor_mass, peaks);
                                spectra.Add(spectrum);
                            }
                        }
                        else//*/
                        {/*
                        if (options.assignChargeStates)
                        {
                            peaks = AssignChargeStatesbkp(peaks, charge, options.productMassTolerance);
                            if (options.deisotope)
                            {
                                peaks = Deisotopebkp(peaks, charge, options.productMassTolerance);
                            }
                        }//*/
                            //peaks = AssignChargeStatesAndDeisotope(peaks, options.maximumAssumedPrecursorChargeState, options.productMassTolerance);

                            double precursor_mass = Numerics.MassFromMZ(precursor_mz, charge);

                            ProductSpectrum spectrum = new ProductSpectrum(scan_number, retention_time, fragmentation_method, precursor_mz, precursor_intensity, charge, precursor_mass, peaks, isolationWindow);
                            spectra.AddMSMS(spectrum);
                            //zones.Add(new Zone(precursor_mz - isolationWindow, precursor_mz + isolationWindow, retention_time));
                        }

                        //if (spectra.Count >= maxNbMSMS)
                        //    i = 10000000;
                    }
                    else //Is an MS
                    {
                        if (loadMS)
                        {
                            double retention_time = spec.scanList.scans[0].cvParam(pwiz.CLI.cv.CVID.MS_scan_start_time).timeInSeconds() / 60.0;

                            pwiz.CLI.msdata.BinaryDataArray mz = spec.getMZArray();
                            pwiz.CLI.msdata.BinaryDataArray intensity = spec.getIntensityArray();

                            if (previousMS1 != null)
                            {
                                previousMS1.ScanDuration = retention_time - previousMS1.RetentionTimeInMin;
                                spectra.MS1s.Add(previousMS1);
                            }
                            previousMS1 = new MS1Spectrum(i, retention_time, intensity.data, mz.data, 1);
                            //Trail.Follow(mz.data, intensity.data, retention_time, ref trails, options);
                            //Trail.RemoveFinished(ref trails, spectra, 1);
                        }
                    }
                    spec.Dispose(); spec = null;
                    Console.Write("\r{0}%   ", ((100 * i) / num_spectra));
                }
                if (previousMS1 != null)
                    spectra.MS1s.Add(previousMS1);
                /*
                //Optimization of Track following parameters
                long nbChargedTracks = 0;
                for(int missingScans = 1; missingScans < 5; missingScans++)
                {
                    for(int centroid = 1; centroid < 5; centroid++)
                    {
                        for(int minPeaks = 1; minPeaks < 7; minPeaks++)
                        {
                            for(double valleyFactor = 0.1; valleyFactor < 4; valleyFactor += 0.3)
                            {                            
                                //weightedMean
                                Tracks tracks = ComputeSpectraTracks(spectra, options, mzMlFilepath, missingScans, centroid, minPeaks, valleyFactor, MaxQuant.CentroidPosition.weightedMean);
                                tracks.Sort(Tracks.AscendingPrecursorMassComparison);
                                long cumulIsotopes = 0;
                                foreach (stTrack track in tracks)
                                    cumulIsotopes += Queries.GetIsotopes(track, options, tracks, sample).Count;
                                if (cumulIsotopes > nbChargedTracks)
                                {
                                    nbChargedTracks = cumulIsotopes;
                                    Console.WriteLine(missingScans + "," + centroid + "," + minPeaks + "," + valleyFactor + ",weightedMean");
                                }
                            
                                //Gaussian
                                tracks = ComputeSpectraTracks(spectra, options, mzMlFilepath, missingScans, centroid, minPeaks, valleyFactor, MaxQuant.CentroidPosition.gaussian);
                                tracks.Sort(Tracks.AscendingPrecursorMassComparison);
                                cumulIsotopes = 0;
                                foreach (stTrack track in tracks)
                                    cumulIsotopes += Queries.GetIsotopes(track, options, tracks, sample).Count;
                                if (cumulIsotopes > nbChargedTracks)
                                {
                                    nbChargedTracks = cumulIsotopes;
                                    Console.WriteLine(missingScans + "," + centroid + "," + minPeaks + "," + valleyFactor + ",Gaussian");
                                }
                            }
                        }
                    }
                }//*/

                if (spectra.MS1s.Count > 0)
                    spectra.tracks = ComputeSpectraTracks(spectra, options, mzMlFilepath, 3, 1, 3, 1.7, MaxQuant.CentroidPosition.weightedMean);
                else
                    spectra.tracks = new Tracks();
                spectra.tracks.Sort(Tracks.AscendingPrecursorMassComparison);
                Console.Write("\r{0}%   ", 100);

                //ContaminantMasses.DisplayContaminants();
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.StackTrace);
                Console.WriteLine(ex.Message);
            }
            return spectra;
        }

        private static Tracks ComputeSpectraTracks(Spectra spectra, DBOptions options, string outputFilename, int missingScan, int centroid, int minPeaks, double valleyFactor, MaxQuant.CentroidPosition centroidMethod)
        {
            //Trail.RemoveFinished(ref trails, spectra, -1);
            double[] centerMassArray;
            float[] centerMassErrorArray;
            float[] intensityArray;
            float[] minTimeArray;
            float[] maxTimeArray;
            long[] filePosArray;
            //TODO Cycle values to optimize missing scans and centroid values

            string file = null;
            if (options.WriteMaxQuantPeakFile)
                file = options.OutputFolder + vsCSV.GetFileName_NoExtension(outputFilename) + "_Peaks.txt";
            MaxQuant.PeakDetection.Detect(  file,
                                            missingScan,//*1-2-3-4-5
                                            centroid,//*1-2-3-4-5-6-7-8-9-10
                                            centroidMethod,//*
                                            false, 0, options.precursorMassTolerance.Value,//TODO ensure its always in ppm
                                            minPeaks,//*1-2-3-4-5-6-7-8-9-10
                                            valleyFactor,//*0.1-0.2-0.3-...-3.0
                                            true,
                                            0,
                                            new Trinity.MaxQuant.RawFileWrapper(spectra),
                                            true,
                                            null, out centerMassArray, out centerMassErrorArray, out intensityArray, out minTimeArray, out maxTimeArray, out filePosArray);

            Tracks tracks = new Tracks();
            for (int i = 0; i < centerMassArray.Length; i++)
                tracks.AddTrack(centerMassArray[i], (minTimeArray[i] + maxTimeArray[i]) * 0.5, minTimeArray[i], maxTimeArray[i], intensityArray[i]);
            return tracks;
        }

        private static GraphML_List<MsMsPeak> AssignChargeStates(GraphML_List<MsMsPeak> peaks, int maxCharge, MassTolerance isotopicMzTolerance)
        {
            GraphML_List<MsMsPeak> new_peaks = new GraphML_List<MsMsPeak>(peaks);

            for (int i = 0; i < peaks.Count - 1; i++)
            {
                double massMax = (peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE) + isotopicMzTolerance;
                int j = i + 1;
                List<int> charges = new List<int>();
                while (j < peaks.Count)
                {
                    if (peaks[j].MZ > massMax)
                        break;

                    for (int c = maxCharge; c >= 1; c--)
                    {
                        if (Math.Abs(Numerics.CalculateMassError(peaks[j].MZ, peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE / (double)c, isotopicMzTolerance.Units)) <= isotopicMzTolerance.Value)
                        {
                            new_peaks.Add(new MsMsPeak(peaks[i].MZ, peaks[i].Intensity, c));
                            charges.Add(c);
                        }
                    }

                    j++;
                }
                if (charges.Count == 0)
                {
                    new_peaks.Add(new MsMsPeak(peaks[i].MZ, peaks[i].Intensity, 0));
                }
            }

            return new_peaks;
        }

        private static GraphML_List<MsMsPeak> AssignChargeStatesbkp(GraphML_List<MsMsPeak> peaks, int maxCharge, MassTolerance isotopicMzTolerance)
        {
            GraphML_List<MsMsPeak> new_peaks = new GraphML_List<MsMsPeak>();

            for(int i = 0; i < peaks.Count - 1; i++)
            {
                int j = i + 1;
                List<int> charges = new List<int>();
                while(j < peaks.Count)
                {
                    if(peaks[j].MZ > (peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE) + isotopicMzTolerance)
                    {
                        break;
                    }

                    for(int c = maxCharge; c >= 1; c--)
                    {
                        if(Math.Abs(Numerics.CalculateMassError(peaks[j].MZ, peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE / (double)c, isotopicMzTolerance.Units)) <= isotopicMzTolerance.Value)
                        {
                            new_peaks.Add(new MsMsPeak(peaks[i].MZ, peaks[i].Intensity, c));
                            charges.Add(c);
                        }
                    }

                    j++;
                }
                if(charges.Count == 0)
                {
                    new_peaks.Add(new MsMsPeak(peaks[i].MZ, peaks[i].Intensity, 0));
                }
            }

            return new_peaks;
        }

        private static GraphML_List<MsMsPeak> AssignChargeStatesAndDeisotope(GraphML_List<MsMsPeak> peaks, int maxCharge, MassTolerance isotopicMzTolerance)
        {
            GraphML_List<MsMsPeak> new_peaks = new GraphML_List<MsMsPeak>(peaks);
            //peaks.Sort(MSPeak.AscendingMzComparison);

            int[] bestIsotopes = new int[4];
            int[] currentIsotopes = new int[4];
            for (int lowMassIndex = 0; lowMassIndex < new_peaks.Count - 1; lowMassIndex++)
            {
                double bestChargeScore = 0;
                int bestCharge = 0;
                bestIsotopes[0] = 0; bestIsotopes[1] = 0; bestIsotopes[2] = 0; bestIsotopes[3] = 0;
                for(int charge = maxCharge; charge > 0; charge--)
                {
                    currentIsotopes[0] = 0; currentIsotopes[1] = 0; currentIsotopes[2] = 0; currentIsotopes[3] = 0;
                    double score = 0;
                    int potentialIsotopeIndex = lowMassIndex + 1;
                    for(int isotope = 1; isotope <= 4; isotope++)
                    {
                        double bestMassError = isotopicMzTolerance.Value;
                        double aim = Numerics.IsotopicMassShift(isotope, charge) + new_peaks[lowMassIndex].MZ;

                        while (potentialIsotopeIndex < new_peaks.Count && new_peaks[potentialIsotopeIndex].MZ < aim + bestMassError)
                        {
                            if (new_peaks[lowMassIndex].Intensity > new_peaks[potentialIsotopeIndex].Intensity)
                            {
                                double massError = Math.Abs(Numerics.CalculateMassError(new_peaks[potentialIsotopeIndex].MZ, aim, isotopicMzTolerance.Units));
                                if (massError < bestMassError)
                                {
                                    bestMassError = massError;
                                    currentIsotopes[isotope-1] = potentialIsotopeIndex;
                                }
                            }
                            potentialIsotopeIndex++;
                        }
                        score += isotopicMzTolerance.Value - bestMassError;
                        if (score == 0)
                            break;;
                    }
                    if (score > bestChargeScore)
                    {
                        bestIsotopes[0] = currentIsotopes[0];
                        bestIsotopes[1] = currentIsotopes[1];
                        bestIsotopes[2] = currentIsotopes[2];
                        bestIsotopes[3] = currentIsotopes[3];
                        bestChargeScore = score;
                        bestCharge = charge;
                    }
                }

                new_peaks[lowMassIndex].Charge = bestCharge;
                for(int i = 3; i >= 0; i--)
                    if (bestIsotopes[i] > 0)
                    {
                        new_peaks[lowMassIndex].Intensity += new_peaks[bestIsotopes[i]].Intensity;
                        new_peaks.RemoveAt(bestIsotopes[i]);
                    }                
            }
            return new_peaks;
        }

        private static GraphML_List<MsMsPeak> Deisotopebkp(GraphML_List<MsMsPeak> peaks, int maxCharge, MassTolerance isotopicMzTolerance)
        {
            GraphML_List<MsMsPeak> new_peaks = new GraphML_List<MsMsPeak>(peaks);
            peaks.Sort(MsMsPeak.AscendingMzComparison);

            for(int lowMassIndex = 0; lowMassIndex < new_peaks.Count - 1; lowMassIndex++)
            {
                if(new_peaks[lowMassIndex].Charge > 0)
                {
                    int toRemove = -1;
                    double bestMassError = isotopicMzTolerance.Value;
                    double aim = Numerics.IsotopicMassShift(1, new_peaks[lowMassIndex].Charge) + new_peaks[lowMassIndex].MZ;

                    int potentialIsotopeIndex = lowMassIndex + 1;
                    while (potentialIsotopeIndex < new_peaks.Count && new_peaks[potentialIsotopeIndex].MZ < aim + bestMassError)
                    {
                        if(new_peaks[lowMassIndex].Intensity > new_peaks[potentialIsotopeIndex].Intensity)
                        {
                            double massError = Math.Abs(Numerics.CalculateMassError(new_peaks[potentialIsotopeIndex].MZ, aim, isotopicMzTolerance.Units));
                            if(massError < bestMassError)
                            {
                                bestMassError = massError;
                                toRemove = potentialIsotopeIndex;
                            }
                        }
                        potentialIsotopeIndex++;
                    }
                    if(toRemove > 0)
                    {
                        new_peaks[lowMassIndex].Intensity += new_peaks[toRemove].Intensity;
                        new_peaks.RemoveAt(toRemove);
                    }
                }
            }
            return new_peaks;
        }

        private static GraphML_List<MsMsPeak> FilterPeaks(GraphML_List<MsMsPeak> peaks, int maximumNumberOfPeaks)
        {
            GraphML_List<MsMsPeak> filtered_peaks = new GraphML_List<MsMsPeak>(peaks);

            if (maximumNumberOfPeaks > 0 && filtered_peaks.Count > maximumNumberOfPeaks)
            {
                filtered_peaks.Sort(MsMsPeak.DescendingIntensityComparison);
                filtered_peaks.RemoveRange(maximumNumberOfPeaks, filtered_peaks.Count - maximumNumberOfPeaks);
            }

            return filtered_peaks;
        }
    }
}