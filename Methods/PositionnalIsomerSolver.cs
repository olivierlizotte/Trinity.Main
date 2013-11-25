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

namespace Trinity.UnitTest
{
    public class CharacterizedPrecursor : PrecursorIon
    {
        public Peptide Peptide;
        public PeptideSpectrumMatches Psms;
        public double NormalizeFactor;
        public List<ProductMatch> Fragments;

        private Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor;        

        public CharacterizedPrecursor(Sample sample, DBOptions dbOptions, Peptide peptide, IEnumerable<Query> queries, double mz, int charge = -1):base(sample, queries, mz, charge)
        {
            this.Peptide = peptide;
            Psms = new PeptideSpectrumMatches();
            foreach(Query query in queries)
                if(query.sample == sample)
                    Psms.Add(new PeptideSpectrumMatch(query, peptide, dbOptions));
            Psms.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
            Fragments = Psms.GetCombinedSpectrum(dbOptions);
            Fragments.Sort(ProductMatch.DescendingWeightComparison);

            DicOfPsmFactor = new Dictionary<PeptideSpectrumMatch, double>();
            foreach (PeptideSpectrumMatch psm in Psms)
                DicOfPsmFactor.Add(psm, psm.ProbabilityScore());
        }

        public List<ProductMatch> GetCombinedMatches(Dictionary<double, int> dicOfCommonFragments, DBOptions dbOptions)
        {
            List<ProductMatch> matches = new List<ProductMatch>(dicOfCommonFragments.Count);

            foreach (double mz in dicOfCommonFragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in Fragments)
                    if (match.theoMz == mz)
                    {
                        matches.Add(new ProductMatch(match));
                        found = true;
                    }

                if (!found)
                {
                    double sumPsmFactor = 0;
                    ProductMatch newMatch = new ProductMatch();
                    newMatch.theoMz = mz;
                    newMatch.weight = 0;
                    newMatch.obsIntensity = 0;
                    newMatch.normalizedIntensity = 0;
                    foreach (PeptideSpectrumMatch psm in Psms)
                    {
                        foreach (MsMsPeak peak in psm.Query.spectrum.Peaks)
                        {
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(peak.MZ, mz, dbOptions.productMassTolerance.Units)) <= dbOptions.productMassTolerance.Value)
                            {
                                newMatch.obsIntensity += peak.Intensity * DicOfPsmFactor[psm];
                                newMatch.normalizedIntensity += (peak.Intensity / (psm.Query.spectrum.PrecursorIntensityPerMilliSecond * psm.Query.spectrum.InjectionTime)) * DicOfPsmFactor[psm];
                                sumPsmFactor += DicOfPsmFactor[psm];
                                newMatch.weight++;
                            }
                        }
                    }
                    if (newMatch.weight > 0)
                    {
                        newMatch.obsIntensity /= sumPsmFactor;
                        newMatch.normalizedIntensity /= sumPsmFactor;
                    }
                    newMatch.weight *= newMatch.normalizedIntensity;
                    matches.Add(newMatch);
                }
            }

            double averageNormedIntensity = 0.0;
            foreach(ProductMatch match in matches)
                averageNormedIntensity += match.normalizedIntensity;

            if (matches.Count > 0)
                averageNormedIntensity /= (double)matches.Count;

            //Keep only most intense fragments (5% of average normalized intensity)
            foreach (ProductMatch pm in matches)
                if (pm.normalizedIntensity < averageNormedIntensity * 0.1)//0.05
                {
                    pm.normalizedIntensity = 0;
                    pm.obsIntensity = 0;
                }
                else
                {
                    pm.normalizedIntensity /= NormalizeFactor;
                    pm.obsIntensity /= NormalizeFactor;
                }
            return matches;
        }

        public void Normalize(double average)
        {
            //average = area * Norm => Norm = average/area
            if(eCurve.Area > 0)
                NormalizeFactor = average / eCurve.Area;
        }
    }

    public class MixedPrecursor : PrecursorIon
    {
        public Dictionary<Sample, ElutionCurve> PeptideRatios;
        public MixedPrecursor(Sample sample, IEnumerable<Query> queries, double mz):base(sample, queries, mz, -1)
        {
        }

        public Dictionary<Sample, double> ComputePeptideRatios(Dictionary<Dictionary<Sample, ElutionCurve>, double> DicOfCurveErrors)
        {
            //!!!! Remove very bad Ratios before going forward
            Remove Bad Ratios
            int nbRun = 3;
            Dictionary<Dictionary<Sample, ElutionCurve>, double> dicOfCorrelations = new Dictionary<Dictionary<Sample,ElutionCurve>,double>(); 
            foreach (Dictionary<Sample, ElutionCurve> dicOfCurve in DicOfCurveErrors.Keys)
                dicOfCorrelations.Add(dicOfCurve, 1.0 / (double) DicOfCurveErrors.Count);

            Dictionary<Sample, double> bestAverage = null;
            while (nbRun > 0)
            {
                nbRun--;

                Dictionary<Sample, double> average = new Dictionary<Sample, double>();
                foreach (Dictionary<Sample, ElutionCurve> dicOfCurve in DicOfCurveErrors.Keys)
                {
                    Dictionary<Sample, double> areas = GetAreas(dicOfCurve);
                    foreach (Sample sample in areas.Keys)
                    {
                        if (!average.ContainsKey(sample))
                            average.Add(sample, 0);
                        average[sample] += areas[sample] * dicOfCorrelations[dicOfCurve];
                    }
                }

                foreach (Dictionary<Sample, ElutionCurve> dicOfCurve in DicOfCurveErrors.Keys)
                {
                    Dictionary<Sample, double> elution = new Dictionary<Sample, double>();
                    foreach (Sample sample in average.Keys)
                        if (dicOfCurve.ContainsKey(sample))
                            elution.Add(sample, dicOfCurve[sample].Area);
                        else
                            elution.Add(sample, 0);
                    dicOfCorrelations[dicOfCurve] = Math.Abs(MathNet.Numerics.Statistics.Correlation.Pearson(average.Values, elution.Values));
                }

                //Normalize correlation factors (sum must equal 1)
                double sumOfCorr = 0.0;
                foreach (Dictionary<Sample, ElutionCurve> dicOfCurve in DicOfCurveErrors.Keys)
                    sumOfCorr += dicOfCorrelations[dicOfCurve];

                foreach (Dictionary<Sample, ElutionCurve> dicOfCurve in DicOfCurveErrors.Keys)
                    dicOfCorrelations[dicOfCurve] /= sumOfCorr;

                bestAverage = average;
            }
            return bestAverage;
        }

        public static Dictionary<Sample, double> GetNormalizedAreas(Dictionary<Sample, ElutionCurve> curves)
        {
            Dictionary<Sample, double> normedCurves = new Dictionary<Trinity.Sample, double>();
            double min = double.MaxValue;
            double max = double.MinValue;
            foreach (ElutionCurve pepCurve in curves.Values)
            {
                if (pepCurve.Area < min)
                    min = pepCurve.Area;

                if (pepCurve.Area > max)
                    max = pepCurve.Area;
            }
            foreach (Sample sample in curves.Keys)
                normedCurves.Add(sample, (curves[sample].Area - min) / (max - min));

            return normedCurves;
        }

        public static Dictionary<Sample, double> GetAreas(Dictionary<Sample, ElutionCurve> curves)
        {
            Dictionary<Sample, double> newCurves = new Dictionary<Trinity.Sample, double>();
            foreach (Sample sample in curves.Keys)
                newCurves.Add(sample, curves[sample].Area);

            return newCurves;
        }
    }

    public class PrecursorIon
    {
        public ElutionCurve eCurve;
        public double MZ;
        public int Charge;
        public List<Query> Queries;
        public Sample Sample;
        public PrecursorIon(Sample sample, IEnumerable<Query> queries, double mz, int charge)
        {
            this.MZ = mz;
            this.Charge = charge;
            this.Queries = new List<Query>();
            foreach (Query query in queries)
                if (query.sample == sample)
                    this.Queries.Add(query);
            this.Queries.Sort(Query.AscendingRetentionTimeComparison);            
            this.eCurve = ElutionCurve.Create(this.Queries);
            this.Sample = sample;
        }
    }

    public class PositionnalIsomerSolver
    {
        private static DBOptions CreateOptions(string fastaFile, string outputFolder, IConSol consol)
        {
            DBOptions dbOptions = new DBOptions(fastaFile, consol);
            dbOptions.precursorMassTolerance = new MassTolerance(8, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(20, MassToleranceUnits.ppm);
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

        public static void Solve(string[] spikedRaws, string[] mixedRaws, string fastaFile, string outputFolder, IConSol conSol)
        {
            DBOptions dbOptions = CreateOptions(fastaFile, outputFolder, conSol);
            Samples SpikedSamples = new Samples(dbOptions);
            for (int i = 0; i < spikedRaws.Length; i++)
                SpikedSamples.Add(new Sample(i + 1, 1, 1, spikedRaws[i], spikedRaws[i], 0, ""));

            //Precompute Spiked peptide identifications
            Result spikedResult = Propheus.Start(dbOptions, SpikedSamples, false, false, true);
            
            Samples MixedSamples = new Samples(dbOptions);
            for (int i = 0; i < mixedRaws.Length; i++)
                MixedSamples.Add(new Sample(i + 1, 1, 1, mixedRaws[i], mixedRaws[i], 0, ""));

            //Precompute Mixed peptide identifications
            Result mixedResult = Propheus.Start(dbOptions, MixedSamples, false, false, true);

            //!!!! Compute spiked peptide profiles and list of fragments
            //!!!! Compute Normalization vectors
            Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> characterizedPeptides = GetSpikedPrecursors(SpikedSamples, spikedResult, dbOptions);

            vsCSVWriter writerCumul = new vsCSVWriter(dbOptions.OutputFolder + "Result.csv");
            string title = "Mixed Sample,Precursor";
            foreach(double precursor in characterizedPeptides.Keys)
                foreach(CharacterizedPrecursor charPrec in characterizedPeptides[precursor].Values)
                    title += "," + charPrec.Peptide.Sequence + " Charge " + charPrec.Charge;
            writerCumul.AddLine(title);

            string curveStr = "Polynomial Curve";
            foreach (double precursor in characterizedPeptides.Keys)
                foreach (CharacterizedPrecursor charPrec in characterizedPeptides[precursor].Values)
                    if(charPrec.eCurve.Area > 0 && charPrec.eCurve.Coefficients.Length == 3)
                        curveStr += "," + charPrec.eCurve.Coefficients[0] + "x^2 + " + charPrec.eCurve.Coefficients[1] + "x" + charPrec.eCurve.Coefficients[2];
            writerCumul.AddLine(curveStr);

            //Get the list of precursors to characterize
            foreach (Sample mixedSample in MixedSamples)
            {
                Dictionary<double, MixedPrecursor> mixedPrecursors = GetMixedPrecursors(mixedSample, mixedResult, dbOptions, characterizedPeptides);

                foreach (double keyMz in mixedPrecursors.Keys)
                {
                    // Compute Max Flow for this precursor
                    Dictionary<Sample, double> ratios = GetRatiosFromMaxFlow(dbOptions, characterizedPeptides, mixedPrecursors[keyMz]);
                    
                    string resultStr = vsCSV.GetFileName(mixedSample.sSDF) + "," + keyMz;
                    foreach (double precursor in characterizedPeptides.Keys)
                        foreach (CharacterizedPrecursor charPrec in characterizedPeptides[precursor].Values)
                            if (ratios.ContainsKey(charPrec.Sample))
                                resultStr += "," + ratios[charPrec.Sample];
                            else
                                resultStr += ",0";
                    writerCumul.AddLine(resultStr);

                    // Export results in a file
                    vsCSVWriter writerRatio = new vsCSVWriter(dbOptions.OutputFolder + vsCSV.GetFileName_NoExtension(mixedSample.sSDF) + "_" + keyMz + "MZ.csv");
                    writerRatio.AddLine("Sequence,Area");
                    writerRatio.AddLine("Total," + mixedPrecursors[keyMz].eCurve.Area);
                    foreach (Sample spikedSample in ratios.Keys)
                        writerRatio.AddLine(characterizedPeptides[keyMz][spikedSample].Peptide.Sequence + "," + ratios[spikedSample]);
                    writerRatio.WriteToFile();
                }
            }
            writerCumul.WriteToFile();
        }

        public static Dictionary<double, List<Query>> GetPrecursors(Result result, Sample sample, DBOptions dbOptions, IEnumerable<double> keys)
        {
            Dictionary<double, List<Query>> DicOfSpectrumMasses = new Dictionary<double, List<Query>>();
            foreach (Query query in result.queries)
            {
                if (query.sample == sample)
                {
                    double foundKey = query.spectrum.PrecursorMZ;
                    
                    foreach (double key in keys)
                        if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(query.spectrum.PrecursorMZ, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                            foundKey = key;                                                

                    if (!DicOfSpectrumMasses.ContainsKey(foundKey))
                    {
                        List<Query> listOfSpectrum = new List<Query>();
                        listOfSpectrum.Add(query);
                        DicOfSpectrumMasses.Add(foundKey, listOfSpectrum);
                    }
                    else
                        DicOfSpectrumMasses[foundKey].Add(query);
                }
            }
            return DicOfSpectrumMasses;
        }

        public static Dictionary<double, MixedPrecursor> GetMixedPrecursors(Sample mixedSample, Result mixedResult, DBOptions dbOptions, Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> charPeptides)
        {
            Dictionary<double, List<Query>> DicOfSpectrumMasses = GetPrecursors(mixedResult, mixedSample, dbOptions, charPeptides.Keys);
            Dictionary<double, MixedPrecursor> DicOfMixedPrecursor = new Dictionary<double,MixedPrecursor>();
            foreach (double key in DicOfSpectrumMasses.Keys)
            {
                if (charPeptides.ContainsKey(key))
                {
                    MixedPrecursor mixedPrecursor = new MixedPrecursor(mixedSample, DicOfSpectrumMasses[key], key);
                    DicOfMixedPrecursor.Add(key, mixedPrecursor);
                }
            }
            return DicOfMixedPrecursor;            
        }

        public static Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> GetSpikedPrecursors(Samples spikedSamples, Result spikedResult, DBOptions dbOptions)
        {
            Dictionary<double, Dictionary<Query, int>> mzKeys = new Dictionary<double, Dictionary<Query, int>>();
            foreach(Query query in spikedResult.queries)
            {
                foreach(PeptideSpectrumMatch psm in query.psms)
                {
                    double mz = Numerics.MZFromMass(psm.Peptide.MonoisotopicMass, query.spectrum.PrecursorCharge);
                    if(!mzKeys.ContainsKey(mz))
                    {
                        bool found = false;
                        foreach(double key in mzKeys.Keys)
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(mz, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                            {
                                mz = key;
                                found = true;
                            }
                        if(!found)
                            mzKeys.Add(mz, new Dictionary<Query, int>());
                    }
                    if(!mzKeys[mz].ContainsKey(query))
                        mzKeys[mz].Add(query, 1);
                }
            }
                
            Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> spikes = new Dictionary<double,Dictionary<Sample,CharacterizedPrecursor>>();
            foreach(Sample spikedSample in spikedSamples)
            {
                Dictionary<double, List<Query>> DicOfSpectrumMasses = GetPrecursors(spikedResult, spikedSample, dbOptions, mzKeys.Keys);                
                foreach(double mzKey in DicOfSpectrumMasses.Keys)
                {
                    if (mzKeys.ContainsKey(mzKey))
                    {
                        //Pick the best PSM for each sample/precursor pair
                        Dictionary<Peptide, double> DicOfProbabilityScores = new Dictionary<Peptide, double>();

                        foreach (Query query in mzKeys[mzKey].Keys)
                            if (query.sample == spikedSample)
                            {
                                foreach (PeptideSpectrumMatch psm in query.psms)
                                    if (!DicOfProbabilityScores.ContainsKey(psm.Peptide))
                                        DicOfProbabilityScores.Add(psm.Peptide, psm.ProbabilityScore());
                                    else
                                        DicOfProbabilityScores[psm.Peptide] += psm.ProbabilityScore();
                            }

                        Peptide bestPeptide = null;
                        double bestScore = double.MinValue;
                        foreach (Peptide keyPep in DicOfProbabilityScores.Keys)
                            if (DicOfProbabilityScores[keyPep] > bestScore)
                            {
                                bestScore = DicOfProbabilityScores[keyPep];
                                bestPeptide = keyPep;
                            }
                        if (bestPeptide != null)
                        {
                            CharacterizedPrecursor cPrec = new CharacterizedPrecursor(spikedSample, dbOptions, bestPeptide, mzKeys[mzKey].Keys, mzKey);
                            if (!spikes.ContainsKey(mzKey))
                                spikes.Add(mzKey, new Dictionary<Sample, CharacterizedPrecursor>());
                            if (!spikes[mzKey].ContainsKey(spikedSample))
                                spikes[mzKey].Add(spikedSample, cPrec);
                            else
                                Console.WriteLine("Twice??");
                        }
                    }
                }//End of foreach mzKey
            }//End of foreach spiked sample

            //Normalize intensities based on average area of each precursor
            foreach (double mzKey in spikes.Keys)
            {
                double cumulArea = 0.0;
                int nbNonZero = 0;
                foreach (CharacterizedPrecursor precursor in spikes[mzKey].Values)
                {
                    if(precursor.eCurve.Area > 0)
                        nbNonZero ++;
                    cumulArea += precursor.eCurve.Area;
                }

                double average = cumulArea / (double)nbNonZero;
                foreach (CharacterizedPrecursor precursor in spikes[mzKey].Values)
                    precursor.Normalize(average);
            }
            return spikes;
        }

        public static Dictionary<double, int> GetCommonFragmentMz(IEnumerable<CharacterizedPrecursor> precursors, int nbFragmentPerPrec)
        {
            Dictionary<double, int> DicOfFragments = new Dictionary<double,int>();
            foreach(CharacterizedPrecursor cPrec in precursors)
            {
                for(int i = 0; i < nbFragmentPerPrec && i < cPrec.Fragments.Count; i++)
                    if(!DicOfFragments.ContainsKey(cPrec.Fragments[i].theoMz))
                        DicOfFragments.Add(cPrec.Fragments[i].theoMz, 1);
                    else
                        DicOfFragments[cPrec.Fragments[i].theoMz] ++;
            }
            return DicOfFragments;
        }
        
        public static Dictionary<Sample, double> GetRatiosFromMaxFlow(DBOptions dbOptions, Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> spikes, MixedPrecursor mixedPrecursor)        
        {
            int minProducts = 4;
            int maxProducts = 14;
            Dictionary<Dictionary<Sample, ElutionCurve>, double> DicOfCurveErrors = new Dictionary<Dictionary<Sample, ElutionCurve>, double>();
            
            for (int nbProductsToKeep = minProducts; nbProductsToKeep <= maxProducts; nbProductsToKeep++)
            {
                Dictionary<Sample, CharacterizedPrecursor> Isomers = new Dictionary<Sample, CharacterizedPrecursor>();
                foreach (double mz in spikes.Keys)
                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(mz, mixedPrecursor.MZ, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                        foreach (Sample sample in spikes[mz].Keys)
                            Isomers.Add(sample, spikes[mz][sample]);

                Dictionary<double, int> dicOfCommonFragments = GetCommonFragmentMz(Isomers.Values, nbProductsToKeep);
                Dictionary<Sample, List<ProductMatch>> matches = new Dictionary<Sample, List<ProductMatch>>();
                foreach (CharacterizedPrecursor cPrec in Isomers.Values)
                    matches.Add(cPrec.Sample, cPrec.GetCombinedMatches(dicOfCommonFragments, dbOptions));

                double cumulUnderflow = 0;
                Dictionary<Sample, ElutionCurve> curves = new Dictionary<Sample,ElutionCurve>();
                foreach (Query query in mixedPrecursor.Queries)
                {
                    double timeInMilliSeconds = query.spectrum.RetentionTimeInMin * 60.0 * 1000.0;
                    double overFlow = 0;
                    double underFlow = 0;
                    double percentError = 0;
                    Dictionary<Sample, double> finalRatios = LaunchMaxFlowFromSpectrum(matches, 1000, query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                            0, query.spectrum.PrecursorIntensityPerMilliSecond * query.spectrum.InjectionTime, ref overFlow, ref underFlow, ref percentError, dbOptions.ConSole);

                    cumulUnderflow += underFlow;
                    if (percentError < 0.5)
                    {
                        //UsedPrecursorArea += localArea;
                        //Dictionary<Sample, double> avgQuantifiedRatios = new Dictionary<Sample, double>();

                        foreach (Sample sRatio in finalRatios.Keys)
                        {
                            if (!curves.ContainsKey(sRatio))
                                curves.Add(sRatio, new ElutionCurve());

                            curves[sRatio].AddPoint(timeInMilliSeconds, finalRatios[sRatio] * query.spectrum.PrecursorIntensityPerMilliSecond);
                        }
                    }
                    else
                        Console.WriteLine("Ignored Spectrum : " + percentError);
                }//End of foreach query

                mixedPrecursor.PeptideRatios = curves;

                DicOfCurveErrors.Add(curves, cumulUnderflow);
                foreach(Sample sRatio in curves.Keys)
                    curves[sRatio].Compute();
                /*
                    if (!values.ContainsKey(sRatio))
                        values.Add(sRatio, curves[sRatio].Area);
                    else
                        values[sRatio] += curves[sRatio].Area;
                }//*/
            }//End of for each nbProduct            

            Dictionary<Sample, double> averagedValues = mixedPrecursor.ComputePeptideRatios(DicOfCurveErrors);
            /*
            foreach (Dictionary<Sample, ElutionCurve> curves in DicOfCurveErrors.Keys)
            {
                foreach (Sample sample in curves.Keys)
                    if (!averagedValues.ContainsKey(sample))
                        averagedValues.Add(sample, curves[sample].Area);
                    else
                        averagedValues[sample] += curves[sample].Area;
            }//*/
            
            

            /*
            //Get the best description of the ratios, based on the multiple runs
            double maxUnderflow = 0.0;
            foreach (double val in DicOfCurveErrors.Values)
                if(val > maxUnderflow)
                    maxUnderflow = val;

            double sumOfUnderflow = 0.0;
            foreach (double val in DicOfCurveErrors.Values)
                sumOfUnderflow += 1.0 - (val / maxUnderflow);

            if (sumOfUnderflow <= 0)
                sumOfUnderflow = 1;

            //Weighted average (the most under flow gets the least weight
            Dictionary<Sample, double> averagedValues = new Dictionary<Sample, double>();
            foreach(Dictionary<Sample, ElutionCurve> curves in DicOfCurveErrors.Keys)
            {
                foreach(Sample sample in curves.Keys)
                    if(!averagedValues.ContainsKey(sample))
                        averagedValues.Add(sample, curves[sample].Area * ((1.0 - (DicOfCurveErrors[curves] / maxUnderflow)) / sumOfUnderflow));
                    else
                        averagedValues[sample] += curves[sample].Area * ((1.0 - (DicOfCurveErrors[curves] / maxUnderflow)) / sumOfUnderflow);
            }//*/

            return averagedValues;
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
                            ConSole.WriteLine("Potential problem with selected fragment masses ");
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
                        if(bestIndexes[index] >= worstFlowRate)
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
            for(int i = 0; i < ratios.Count; i++)
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
                if(virtualFrag[(float)match.theoMz] > 0)
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

        private static Dictionary<Sample, double> LaunchMaxFlowFromSpectrum(Dictionary<Sample, List<ProductMatch>> ratiosToFit,  
                                            int precision, GraphML_List<MsMsPeak> capacity, MassTolerance tolerance, 
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
                double intensity = peak.Intensity;
                expandedCapacity.Add(new MsMsPeak(peak.MZ, intensity * precision, peak.Charge));
                sumOfProducts += peak.Intensity;
            }
            List<List<ProductMatch>> tmpRatiosToFit = new List<List<ProductMatch>>();
            foreach (List<ProductMatch> list in ratiosToFit.Values)
            {
                List<ProductMatch> pms = new List<ProductMatch>();
                foreach (ProductMatch pm in list)
                {
                    ProductMatch newPm = new ProductMatch(pm);
                    newPm.obsIntensity = newPm.normalizedIntensity * PrecursorIntensityInCTrap;
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

            List<double> result = null;
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
            Dictionary<Sample, double> resultPerSample = new Dictionary<Sample, double>();
            int i = 0;
            foreach (Sample key in ratiosToFit.Keys)
            {
                resultPerSample.Add(key, result[i]);
                i++;
            }
            return resultPerSample;
        }

        private static List<double> GetResultList(List<double> solution, int precision, double underFlow, double sumOfIntensities)
        {
            List<double> rez = new List<double>();
            double sumVal = 0.0;
            foreach(double val in solution)
                sumVal += val;
            
            sumVal += (underFlow / sumOfIntensities) * precision;

            foreach (double val in solution)
                rez.Add(val / sumVal);

            return rez;
        }
    }
}
