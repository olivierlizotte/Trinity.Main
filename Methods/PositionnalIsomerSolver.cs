/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Proteomics.Utilities;
using Trinity.Structures.PositionnalIsomer;

namespace Trinity.Methods
{
    public class PositionnalIsomerSolver
    {
        public int nbMinFragments = 5;
        public int nbMaxFragments = 5;
        public double precTolPpm = 8;
        public double prodTolPpm = 20;
        public long precision     = 1000;
        private DBOptions dbOptions;
        
        private Samples SpikedSamples;
        private Result  SpikedResult;
        private Samples MixedSamples;
        private Result  mixedResult;
        public Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> characterizedPeptides;
        public Dictionary<Sample, Dictionary<double, MixedPrecursor>>         mixedPrecursors;

        private DBOptions CreateOptions(string fastaFile, string outputFolder, IConSol consol)
        {
            DBOptions dbOptions = new DBOptions(fastaFile, consol);
            dbOptions.precursorMassTolerance = new MassTolerance(precTolPpm, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(prodTolPpm, MassToleranceUnits.ppm);
            dbOptions.MaximumPeptideMass = 200000;
            dbOptions.OutputFolder = outputFolder;

            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            dbOptions.DigestionEnzyme = proteases["no enzyme"];
            dbOptions.NoEnzymeSearch = false;
            dbOptions.DecoyFusion = false;
            dbOptions.MaximumNumberOfFragmentsPerSpectrum = 400;
            dbOptions.ToleratedMissedCleavages = 200;
            dbOptions.MinimumPeptideLength = 5;
            dbOptions.MaximumPeptideLength = 300;

            GraphML_List<Modification> fixMods = new GraphML_List<Modification>();
            dbOptions.fixedModifications = fixMods;

            GraphML_List<Modification> varMods = new GraphML_List<Modification>();
            foreach (string strMod in ModificationDictionary.Instance.Keys)
                varMods.Add(ModificationDictionary.Instance[strMod]);

            dbOptions.maximumVariableModificationIsoforms = 1024;
            dbOptions.variableModifications = varMods;

            dbOptions.addFragmentLoss = false;
            dbOptions.addFragmentMods = false;
            dbOptions.fragments = new Fragments();

            dbOptions.fragments.Add(new FragmentA());
            dbOptions.fragments.Add(new FragmentB());
            dbOptions.fragments.Add(new FragmentC());
            dbOptions.fragments.Add(new FragmentX());
            dbOptions.fragments.Add(new FragmentY());
            dbOptions.fragments.Add(new FragmentZ());

            dbOptions.SaveMS1Peaks = true;
            dbOptions.SaveMSMSPeaks = true;
            dbOptions.LoadSpectraIfFound = true;

            dbOptions.NbPSMToKeep = 100;
            return dbOptions;
        }

        public void Solve(string[] spikedRaws, string[] mixedRaws, string fastaFile, string outputFolder, IConSol conSol)
        {
            dbOptions = CreateOptions(fastaFile, outputFolder, conSol);
            SpikedSamples = new Samples(dbOptions);
            for (int i = 0; i < spikedRaws.Length; i++)
                SpikedSamples.Add(new Sample(i + 1, 1, 1, spikedRaws[i], spikedRaws[i], 0, ""));

            //Precompute Spiked peptide identifications
            SpikedResult = Propheus.Start(dbOptions, SpikedSamples, false, false, true);
            
            MixedSamples = new Samples(dbOptions);
            for (int i = 0; i < mixedRaws.Length; i++)
                MixedSamples.Add(new Sample(i + 1, 1, 1, mixedRaws[i], mixedRaws[i], 0, ""));

            //Precompute Mixed peptide identifications
            mixedResult = Propheus.Start(dbOptions, MixedSamples, false, false, true);

            //Compute all usable spiked peptides
            characterizedPeptides = CharacterizedPrecursor.GetSpikedPrecursors(SpikedSamples, SpikedResult, dbOptions, nbMinFragments, nbMaxFragments, precision);
            ExportSpikedSampleResult(characterizedPeptides, dbOptions);

            vsCSVWriter writerCumul = new vsCSVWriter(dbOptions.OutputFolder + @"Combined\Result.csv");
            string titleCombined = "Mixed Sample,Precursor";
            foreach(double precursor in characterizedPeptides.Keys)
                foreach(CharacterizedPrecursor charPrec in characterizedPeptides[precursor].Values)
                    titleCombined += "," + charPrec.Peptide.Sequence + " Charge " + charPrec.Charge;
            writerCumul.AddLine(titleCombined);

            string curveStr = "Polynomial Curve";
            foreach (double precursor in characterizedPeptides.Keys)
                foreach (CharacterizedPrecursor charPrec in characterizedPeptides[precursor].Values)
                    if (charPrec.eCurve.Coefficients != null && charPrec.eCurve.Coefficients.Length == 3)
                        curveStr += "," + charPrec.eCurve.Coefficients[0] + "x^2 + " + charPrec.eCurve.Coefficients[1] + "x" + charPrec.eCurve.Coefficients[2];
                    else
                        curveStr += ",NA";
            writerCumul.AddLine(curveStr);

            mixedPrecursors = new Dictionary<Sample, Dictionary<double, MixedPrecursor>>();

            foreach (Sample mixedSample in MixedSamples) 
                mixedPrecursors.Add(mixedSample, MixedPrecursor.GetMixedPrecursors(mixedSample, mixedResult, dbOptions, characterizedPeptides));

            foreach(double keyMz in characterizedPeptides.Keys)
            {
                //Get the list of precursors to characterize
                foreach (Sample mixedSample in MixedSamples)
                {
                    if (mixedPrecursors[mixedSample].ContainsKey(keyMz))
                    {
                        // Compute Max Flow for this precursor
                        Dictionary<CharacterizedPrecursor, ElutionCurve> ratios = GetRatiosFromMaxFlow(characterizedPeptides, mixedPrecursors[mixedSample][keyMz]);

                        string resultStr = vsCSV.GetFileName(mixedSample.sSDF) + "," + keyMz;
                        foreach (double precursor in characterizedPeptides.Keys)
                            foreach (CharacterizedPrecursor charPrec in characterizedPeptides[precursor].Values)
                                if (ratios.ContainsKey(charPrec))
                                    resultStr += "," + ratios[charPrec].Area;
                                else
                                    resultStr += ",0";
                        writerCumul.AddLine(resultStr);

                        ExportMixedSampleResult(ratios, mixedSample, mixedPrecursors, keyMz, dbOptions);
                    }
                }
            }
            writerCumul.WriteToFile();
        }

        private static void ExportMixedSampleResult(Dictionary<CharacterizedPrecursor, ElutionCurve> ratios, Sample mixedSample, Dictionary<Sample, Dictionary<double, MixedPrecursor>> mixedPrecursors, double keyMz, DBOptions dbOptions)
        {
            // Export results in a file
            vsCSVWriter writerRatio = new vsCSVWriter(dbOptions.OutputFolder + @"Individual\" + vsCSV.GetFileName_NoExtension(mixedSample.sSDF) + "_" + keyMz + "MZ.csv");
            string titleIndividual = "Scan time,Total Area";
            foreach (CharacterizedPrecursor charPep in ratios.Keys)
                titleIndividual += "," + charPep.Peptide.Sequence;

            string line = "Total," + mixedPrecursors[mixedSample][keyMz].eCurve.Area;
            foreach (CharacterizedPrecursor charPep in ratios.Keys)
                line += "," + ratios[charPep].Area;
            writerRatio.AddLine(line);

            for (int i = 0; i < mixedPrecursors[mixedSample][keyMz].eCurve.intensityCount.Count; i++)
            {
                line = mixedPrecursors[mixedSample][keyMz].eCurve.time[i] / (1000.0 * 60.0) + "," + mixedPrecursors[mixedSample][keyMz].eCurve.intensityCount[i];
                foreach (CharacterizedPrecursor charPep in ratios.Keys)
                    line += "," + ratios[charPep].InterpolateIntensity(mixedPrecursors[mixedSample][keyMz].eCurve.time[i]);
                writerRatio.AddLine(line);
            }
            writerRatio.WriteToFile();
        }

        private static void ExportSpikedSampleResult(Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> characterizedPeptides, DBOptions dbOptions)
        {
            foreach (double keyMz in characterizedPeptides.Keys)
            {
                foreach (Sample sample in characterizedPeptides[keyMz].Keys)
                {
                    vsCSVWriter writerRatio = new vsCSVWriter(dbOptions.OutputFolder + @"Individual\" + vsCSV.GetFileName_NoExtension(sample.sSDF) + "_" + keyMz + "MZ.csv");
                    string titleIndividual = "Scan time,Precursor Intensity,Intensity Per Millisecond";
                    foreach (ProductMatch pm in characterizedPeptides[keyMz][sample].AllFragments)
                        titleIndividual += "," + pm.Fragment + pm.fragmentPos + "^" + pm.charge;
                    writerRatio.AddLine(titleIndividual);

                    foreach (Query query in characterizedPeptides[keyMz][sample].Queries)
                    {
                        string line = query.spectrum.RetentionTimeInMin + "," + query.spectrum.PrecursorIntensity + "," + query.spectrum.PrecursorIntensityPerMilliSecond;
                        foreach (ProductMatch pm in characterizedPeptides[keyMz][sample].AllFragments)
                        {
                            double intensity = 0.0;
                            foreach (ProductMatch pmSpec in query.psms[0].AllProductMatches)
                                if (pmSpec.charge == pm.charge && pmSpec.Fragment == pm.Fragment && pmSpec.fragmentPos == pm.fragmentPos)
                                    intensity = pmSpec.obsIntensity;
                            line += "," + intensity;
                        }
                        writerRatio.AddLine(line);
                    }
                    writerRatio.WriteToFile();
                }
            }
        }

        private Dictionary<CharacterizedPrecursor, ElutionCurve> GetRatiosFromMaxFlow(Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> spikes, MixedPrecursor mixedPrecursor)        
        {
            Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double> DicOfCurveErrors = new Dictionary<Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>, double>();
            
            for (int nbProductsToKeep = nbMinFragments; nbProductsToKeep <= nbMaxFragments; nbProductsToKeep++)
            {
                bool validProducts = true;
                int nbIgnoredSpectrum = 0;
                List<CharacterizedPrecursor> Isomers = new List<CharacterizedPrecursor>();
                foreach (double mz in spikes.Keys)
                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(mz, mixedPrecursor.MZ, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                        foreach (Sample sample in spikes[mz].Keys)
                            if (!spikes[mz][sample].Fragments.ContainsKey(nbProductsToKeep))
                                validProducts = false;
                            else
                                Isomers.Add(spikes[mz][sample]);
                if (validProducts)
                {
                    double cumulError = 0;
                    Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> curves = new Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>();
                    foreach (Query query in mixedPrecursor.Queries)
                    {
                        double timeInMilliSeconds = query.spectrum.RetentionTimeInMin * 60.0 * 1000.0;
                        double overFlow = 0;
                        double underFlow = 0;
                        double percentError = 0;
                        Dictionary<CharacterizedPrecursor, MaxFlowResult> finalRatios = LaunchMaxFlowFromSpectrum(Isomers, nbProductsToKeep, 1000, query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                                0, query.spectrum.PrecursorIntensityPerMilliSecond * query.spectrum.InjectionTime, ref overFlow, ref underFlow, ref percentError, dbOptions.ConSole);

                        cumulError += percentError;
                        if (percentError < 0.5)
                        {
                            foreach (CharacterizedPrecursor cPep in finalRatios.Keys)
                            {
                                if (!curves.ContainsKey(cPep))
                                    curves.Add(cPep, new MaxFlowElutionCurve(nbProductsToKeep));

                                //curves[cPep].AddPoint(timeInMilliSeconds, finalRatios[cPep].Ratio * query.spectrum.PrecursorIntensityPerMilliSecond);
                                curves[cPep].AddPoint(timeInMilliSeconds, finalRatios[cPep].Ratio * mixedPrecursor.eCurve.InterpolateIntensity(timeInMilliSeconds));
                            }
                        }
                        else
                            nbIgnoredSpectrum++;

                        if (nbIgnoredSpectrum * 2 > mixedPrecursor.Queries.Count)
                            break;
                    }//End of foreach query

                    if (nbIgnoredSpectrum * 2 < mixedPrecursor.Queries.Count)
                    {
                        if (nbIgnoredSpectrum > 0)
                            Console.WriteLine("Ignored Spectrum : " + nbIgnoredSpectrum + " / " + mixedPrecursor.Queries.Count);

                        foreach (CharacterizedPrecursor cPep in curves.Keys)
                            curves[cPep].Compute();

                        Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve> curvesToKeep = new Dictionary<CharacterizedPrecursor, MaxFlowElutionCurve>();
                        foreach (CharacterizedPrecursor cPep in curves.Keys)
                            if (curves[cPep].Area > 0)
                                curvesToKeep.Add(cPep, curves[cPep]);

                        if (curvesToKeep.Count > 0)
                            DicOfCurveErrors.Add(curvesToKeep, cumulError);
                    }
                }
            }//End of for each nbProduct            

            Dictionary<CharacterizedPrecursor, ElutionCurve> averagedValues = mixedPrecursor.ComputePeptideRatios(DicOfCurveErrors);
            return averagedValues;
        }
        
        public class MaxFlowResult
        {
            public double Ratio;
            public double NbFitTimes;
            public MaxFlowResult(double ratio, double nbTimes)
            {
                this.Ratio = ratio;
                this.NbFitTimes = nbTimes;
            }
        }

        public static Dictionary<CharacterizedPrecursor, MaxFlowResult> LaunchMaxFlowFromSpectrum(IEnumerable<CharacterizedPrecursor> ratiosToFit, int nbProductsToKeep, 
                                            long precision, IEnumerable<MsMsPeak> capacity, MassTolerance tolerance, 
                                            int returnType,//0 for max flow, 1 for best flow, 2 for average
                                            double PrecursorIntensityInCTrap,
                                            ref double overFlow, ref double underFlow, ref double errorInPercent, IConSol ConSole)
        {
            List<List<double>> solutions = new List<List<double>>();
            List<long> average = new List<long>();
            List<MsMsPeak> expandedCapacity = new List<MsMsPeak>();
            double sumOfProducts = 0;
            foreach (MsMsPeak peak in capacity)
            {
                double intensityNormed = peak.Intensity / PrecursorIntensityInCTrap;
                expandedCapacity.Add(new MsMsPeak(peak.MZ, intensityNormed * precision, peak.Charge));
                sumOfProducts += peak.Intensity;
            }
            List<List<ProductMatch>> tmpRatiosToFit = new List<List<ProductMatch>>();
            //foreach (List<ProductMatch> list in ratiosToFit.Values)
            foreach(CharacterizedPrecursor prec in ratiosToFit)
            {
                List<ProductMatch> pms = new List<ProductMatch>();
                foreach (ProductMatch pm in prec.Fragments[nbProductsToKeep])
                {
                    ProductMatch newPm = new ProductMatch(pm);
                    newPm.obsIntensity = newPm.normalizedIntensity;// *PrecursorIntensityInCTrap;
                    pms.Add(newPm);
                }
                tmpRatiosToFit.Add(pms);
            }

            double error = ComputeMaxFlow(tmpRatiosToFit, expandedCapacity, tolerance, ref solutions, ref errorInPercent, ref average, ConSole);

            double sumOfIntensities = 0;
            foreach (MsMsPeak peak in expandedCapacity)
                sumOfIntensities += peak.Intensity;

            overFlow = 0;
            underFlow = error;

            List<MaxFlowResult> result = null;
            switch (returnType)
            {
                case 0:
                    result = GetResultList(solutions[0], precision, underFlow, sumOfIntensities);
                    break;
                case 1:
                    result = GetResultList(solutions[1], precision, underFlow, sumOfIntensities);
                    break;
                case 2:
                    List<double> tmpAverage = new List<double>();                    
                    foreach (double val in average)
                        tmpAverage.Add(val);
                    result = GetResultList(tmpAverage, precision, underFlow, sumOfIntensities);
                    break;
            }
            Dictionary<CharacterizedPrecursor, MaxFlowResult> resultPerSample = new Dictionary<CharacterizedPrecursor, MaxFlowResult>();
            int i = 0;
            foreach (CharacterizedPrecursor key in ratiosToFit)
            {
                resultPerSample.Add(key, result[i]);
                i++;
            }
            return resultPerSample;
        }

        private static List<MaxFlowResult> GetResultList(List<double> solution, long precision, double underFlow, double sumOfIntensities)
        {
            List<MaxFlowResult> rez = new List<MaxFlowResult>();
            double sumVal = 0.0;
            foreach(double val in solution)
                sumVal += val;
            
            sumVal += (underFlow / sumOfIntensities) * precision;

            foreach (double val in solution)
                rez.Add(new MaxFlowResult( val / sumVal, val));

            return rez;
        }

        private static double ComputeMaxFlow(List<List<ProductMatch>> spikedMatches,
                                    List<MsMsPeak> mixedSpectrum, MassTolerance tolerance,
                                ref List<List<double>> optimalSolutions,
                                ref double percentError,
                                ref List<long> average, IConSol ConSole)
        {
            //Create dictionnary of usefull peaks
            Dictionary<float, double> mixedFragDic = new Dictionary<float, double>();
            foreach (List<ProductMatch> fragmentRatio in spikedMatches)
            {
                foreach (ProductMatch match in fragmentRatio)
                {
                    if (!mixedFragDic.ContainsKey((float)match.theoMz))
                    {
                        float closest = -1;
                        foreach (float key in mixedFragDic.Keys)
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(match.theoMz, key, tolerance.Units)) <= tolerance.Value)
                                closest = key;
                        if (closest > 0)
                        {
                            //ConSole.WriteLine("Potential problem with selected fragment masses ");
                            match.theoMz = closest;
                        }
                        else
                            mixedFragDic.Add((float)match.theoMz, 0);
                    }
                }
            }

            //Fill dictionnary with Intensities
            List<float> keys = new List<float>(mixedFragDic.Keys);
            foreach (MsMsPeak peak in mixedSpectrum)
                foreach (float mz in keys)
                {
                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(peak.MZ, mz, tolerance.Units)) <= tolerance.Value)
                        mixedFragDic[mz] += peak.Intensity;
                }

            List<long> localFlows = new List<long>();
            foreach (List<ProductMatch> fragmentRatio in spikedMatches)
                localFlows.Add(FindLocalMaximumFlow(fragmentRatio, mixedFragDic));

            Dictionary<float, double> virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
            double overError = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);
            double underError = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
            double[] bestIndexes = new double[spikedMatches.Count];

            int iterSize = 1;
            double bestOverallError = double.MaxValue;
            List<long> bestLocalFlows = new List<long>();
            Random rnd = new Random();
            while (overError >= 1 && iterSize < 10000)//anything less than 1 is an acceptable solution
            {
                for (int index = 0; index < bestIndexes.Length; index++)
                    bestIndexes[index] = -1;

                for (int i = 0; i < spikedMatches.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= iterSize;

                        virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
                        double tmpErrorMinus = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
                        double tmpErrorPlus = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);

                        double tmpFlowRate = Math.Abs(overError - tmpErrorPlus);
                        double underDiff = Math.Abs(underError - tmpErrorMinus);
                        if (underDiff >= 1)
                            tmpFlowRate /= underDiff;
                        bestIndexes[i] = tmpFlowRate;

                        localFlows[i] += iterSize;
                    }
                }

                //Pick pseudo randomly best index
                double worstFlowRate = 0.0;
                for (int index = 0; index < bestIndexes.Length; index++)
                    if (bestIndexes[index] > worstFlowRate)
                    {
                        worstFlowRate = bestIndexes[index];
                    }

                if (worstFlowRate > 0)
                {
                    int nbMatching = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if (bestIndexes[index] >= worstFlowRate)
                            nbMatching++;

                    int iterChoice = rnd.Next(0, nbMatching - 1);
                    int iterNb = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if (bestIndexes[index] >= worstFlowRate)
                        {
                            if (iterChoice == iterNb)
                                localFlows[index] -= iterSize;
                            iterNb++;
                        }
                    iterSize = 1;
                }
                else
                    iterSize++;

                virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
                overError = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);
                underError = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
                if (overError + underError < bestOverallError)
                {
                    bestLocalFlows = new List<long>(localFlows);
                    bestOverallError = overError + underError;
                }
            }//End of while overflow > 1
            optimalSolutions.Clear();

            List<double> newList = new List<double>();
            foreach (long localFlow in localFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            newList = new List<double>();
            foreach (long localFlow in bestLocalFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            //Compute average
            if (bestOverallError < double.MaxValue)
            {
                average.Clear();
                for (int i = 0; i < optimalSolutions[0].Count; i++)
                {
                    double sum = 0.0;
                    foreach (List<double> solution in optimalSolutions)
                        sum += solution[i];
                    double avg = sum / (double)optimalSolutions.Count;
                    average.Add((long)avg);
                }
            }

            //Compute expected error in percentage
            double sumOfIntensities = 0;
            foreach (double val in mixedFragDic.Values)
                sumOfIntensities += val;
            percentError = underError / sumOfIntensities;

            virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);

            return MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
        }

        private static Dictionary<float, double> BuildVirtualSpectrum(List<List<ProductMatch>> fragRatios, List<long> ratios, Dictionary<float, double> fragments)
        {
            Dictionary<float, double> virtualFrag = new Dictionary<float, double>(fragments.Count);
            for (int i = 0; i < ratios.Count; i++)
            {
                foreach (ProductMatch match in fragRatios[i])
                {
                    if (!virtualFrag.ContainsKey((float)match.theoMz))
                        virtualFrag.Add((float)match.theoMz, match.obsIntensity * ratios[i]);
                    else
                        virtualFrag[(float)match.theoMz] += match.obsIntensity * ratios[i];
                }
            }
            return virtualFrag;
        }

        private static long FindLocalMaximumFlow(List<ProductMatch> spikedMatches, Dictionary<float, double> fragments)
        {
            Dictionary<float, double> virtualFrag = new Dictionary<float, double>(fragments.Count);

            bool letsGo = false;
            double maxCumul = 0;
            foreach (ProductMatch match in spikedMatches)
            {
                if (!virtualFrag.ContainsKey((float)match.theoMz))
                    virtualFrag.Add((float)match.theoMz, match.obsIntensity);
                else
                    virtualFrag[(float)match.theoMz] += match.obsIntensity;
                if (virtualFrag[(float)match.theoMz] > 0)
                    letsGo = true;
                maxCumul += match.obsIntensity;
            }
            maxCumul *= 10000000;
            long nbCumul = 0;
            if (letsGo)
            {
                nbCumul = 1;
                while (nbCumul < maxCumul)
                {
                    foreach (float key in fragments.Keys)
                    {
                        if (virtualFrag.ContainsKey(key))
                        {
                            if (virtualFrag[key] * nbCumul > fragments[key] + 1)
                                return nbCumul;
                        }
                    }
                    nbCumul++;
                }
            }
            return nbCumul;
        }
    }
}
