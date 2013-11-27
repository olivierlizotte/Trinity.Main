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
        public double PrecursorLossNormalizeFactor = 1.0;
        public double FragmentLossNormalizeFactor = 1.0;
        public List<ProductMatch> Fragments;

        private Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor;        

        public CharacterizedPrecursor(Sample sample, DBOptions dbOptions, Peptide peptide, IEnumerable<Query> queries, double mz):base(sample, queries, mz, -1)
        {
            this.Peptide = peptide;
            Psms = new PeptideSpectrumMatches();

            foreach(Query query in queries)
                if (query.sample == sample)
                {
                    Psms.Add(new PeptideSpectrumMatch(query, peptide, dbOptions));
                    Charge = query.precursor.Charge;
                }
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
                //else
                //{
                    //pm.normalizedIntensity /= NormalizeFactor;
                    //pm.obsIntensity /= NormalizeFactor;
                //}
            return matches;
        }

        public double NormalizePrecursor(IEnumerable<CharacterizedPrecursor> allCorrespondingPrec)
        {       
            double average = eCurve.Area;

            //average = area * Norm => Norm = average/area
            if (eCurve.Area > 0)
            {
                //Normalize matches based on precursor intensity differences
                double cumulArea = 0.0;
                int nbNonZero = 0;
                foreach (CharacterizedPrecursor precursor in allCorrespondingPrec)
                {
                    if (precursor.eCurve.Area > 0)
                    {
                        nbNonZero++;
                        cumulArea += precursor.eCurve.Area;
                    }
                }
                average = cumulArea / (double)nbNonZero;
                PrecursorLossNormalizeFactor = Math.Log(average, 2) / Math.Log(this.eCurve.Area, 2);
                //PrecursorLossNormalizeFactor = average / this.eCurve.Area;
            }
            else
                PrecursorLossNormalizeFactor = 1;

            if (PrecursorLossNormalizeFactor > 4) PrecursorLossNormalizeFactor = 4;
            if (PrecursorLossNormalizeFactor < 0.25) PrecursorLossNormalizeFactor = 0.25;
            
            return average;
        }

        public double GetAverageNbTimesFit(DBOptions dbOptions, out double area)
        {        
            long nbRatios = 0;
            double averageTimes = 0;

            int nbAreas = 0;
            area = 0.0;
            for (int nbProductsToKeep = 4; nbProductsToKeep <= 14; nbProductsToKeep++)
            {
                int nbIgnoredSpectrum = 0;
                List<CharacterizedPrecursor> Isomers = new List<CharacterizedPrecursor>();
                Isomers.Add(this);

                Dictionary<double, int> dicOfCommonFragments = PositionnalIsomerSolver.GetCommonFragmentMz(Isomers, nbProductsToKeep);
                Dictionary<CharacterizedPrecursor, List<ProductMatch>> matches = new Dictionary<CharacterizedPrecursor, List<ProductMatch>>();
                foreach (CharacterizedPrecursor cPrec in Isomers)
                    matches.Add(cPrec, cPrec.GetCombinedMatches(dicOfCommonFragments, dbOptions));

                double tmpAverageTimes = 0;
                int tmpNbRatios = 0;
                ElutionCurve curve = new ElutionCurve();
                foreach (Query query in this.Queries)
                {
                    double overFlow = 0;
                    double underFlow = 0;
                    double percentError = 0;
                    Dictionary<CharacterizedPrecursor, PositionnalIsomerSolver.MaxFlowResult> finalRatios = PositionnalIsomerSolver.LaunchMaxFlowFromSpectrum(matches, 1000, query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                            0, query.spectrum.PrecursorIntensityPerMilliSecond * query.spectrum.InjectionTime, ref overFlow, ref underFlow, ref percentError, dbOptions.ConSole);
                    curve.AddPoint(query.spectrum.RetentionTimeInMin * 1000.0 * 60.0, finalRatios[this].Ratio * query.spectrum.PrecursorIntensityPerMilliSecond);
                    if (percentError < 0.5)
                    {
                        tmpAverageTimes += finalRatios[this].NbFitTimes;
                        tmpNbRatios++;
                    }
                }
                if (curve.intensityCount.Count > this.Queries.Count * 0.5)
                {
                    curve.Compute();
                    if (curve.Area > 0)
                    {
                        averageTimes += tmpAverageTimes;
                        nbRatios += tmpNbRatios;
                        area += curve.Area;
                        nbAreas++;
                    }
                }
            }
            area /= (double)nbAreas;
            return averageTimes / (double)nbRatios;
        }

        public void NormalizeFragments(IEnumerable<CharacterizedPrecursor> allCorrespondingPrec, DBOptions dbOptions)
        {
            FragmentLossNormalizeFactor = 1.0;
            if(eCurve.Area > 0)
            {              
                double area = 0.0;
                int nbIter = 3;
                while (nbIter > 0)
                {
                    nbIter--;
                    
                    double averageNbTimes = GetAverageNbTimesFit(dbOptions, out area);
                    FragmentLossNormalizeFactor = 1000.0 / averageNbTimes;
                    foreach (ProductMatch pm in this.Fragments)
                    {
                        pm.normalizedIntensity /= FragmentLossNormalizeFactor;
                        pm.obsIntensity /= FragmentLossNormalizeFactor;
                    }
                }

                double average = NormalizePrecursor(allCorrespondingPrec);
                //PrecursorLossNormalizeFactor = Math.Log(average, 2) / Math.Log(area, 2);


                //PrecursorLossNormalizeFactor = average / area;

                PrecursorLossNormalizeFactor = 1.0;

                    /*
                        //PrecursorLossNormalizeFactor = average / rez[this];
                        //FragmentLossNormalizeFactor = average / rez[this];
                        FragmentLossNormalizeFactor = Math.Log(average, 2) / Math.Log(rez[this], 2);
                        foreach (ProductMatch pm in this.Fragments)
                        {
                            pm.normalizedIntensity *= FragmentLossNormalizeFactor;
                            pm.obsIntensity *= FragmentLossNormalizeFactor;
                        }
                    }//*/
            }
        }//*/
        /*
        public void NormalizeFragments(List<CharacterizedPrecursor> allCorrespondingPrec, Dictionary<CharacterizedPrecursor, List<ProductMatch>> matches)
        {
            //average = area * Norm => Norm = average/area
            if (eCurve.Area > 0)
            {
                //Normalize matches based on normalized fragment intensities
                double cumulFragments = 0.0;
                int nbNonZero = 0;
                foreach (CharacterizedPrecursor precursor in allCorrespondingPrec)
                {
                    if (precursor.eCurve.Area > 0)
                    {
                        foreach (ProductMatch match in matches[precursor])
                            cumulFragments += match.normalizedIntensity;
                        nbNonZero++;
                    }
                }

                double averageFragment = cumulFragments / (double)nbNonZero;
                //NormalizeFactor = average / eCurve.Area;
                //NormalizeFactor = Math.Pow(average, 2) / Math.Pow(eCurve.Area, 2);

                double totalFrag = 0.0;
                foreach (ProductMatch match in matches[this])
                    totalFrag += match.normalizedIntensity;
                FragmentLossNormalizeFactor = averageFragment / totalFrag;
                //NormalizeFactor *= Math.Log(averageFragmentTotal, 2) / Math.Log(totalFrag, 2);

                //double lossOfPrecIntensity = 1;// Math.Log(averagePrecursorArea, 2) / Math.Log(PrecursorCurves[foundKey][sample].Area, 2);
                //double lossofFragmentIntensity = 1;// Math.Pow(avgRatioSum, 2) / Math.Pow(listOfsumOfProducts[sample], 2);

                if (FragmentLossNormalizeFactor > 4) FragmentLossNormalizeFactor = 4;
                if (FragmentLossNormalizeFactor < 0.25) FragmentLossNormalizeFactor = 0.25;

                //foreach (ProductMatch pm in matches[this])
                //{
                //    pm.normalizedIntensity /= NormalizeFactor;
                //    pm.obsIntensity /= NormalizeFactor;
                //}
            }
        }//*/
    }

    public class MixedPrecursor : PrecursorIon
    {
        public Dictionary<Sample, ElutionCurve> PeptideRatios;
        public MixedPrecursor(Sample sample, IEnumerable<Query> queries, double mz):base(sample, queries, mz, -1)
        {
        }

        public Dictionary<CharacterizedPrecursor, double> ComputePeptideRatios(Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double> dicOfCurveErrorsP)
        {
            Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double> dicOfCorrelations = new Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double>();
            foreach (Dictionary<CharacterizedPrecursor, ElutionCurve> dicOfCurve in dicOfCurveErrorsP.Keys)
                dicOfCorrelations.Add(dicOfCurve, 1.0 / (double)dicOfCurveErrorsP.Count);

            int nbRun = 1;

            Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double> lastDicOfCurves = dicOfCurveErrorsP;
            while (nbRun > 0)
            {
                nbRun--;

                //Purge worst curves
                Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double> dicOfCurves = new Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double>();
                double median = MathNet.Numerics.Statistics.Statistics.Median(lastDicOfCurves.Values);
                double maxMed = median + 0.5 * MathNet.Numerics.Statistics.Statistics.Variance(lastDicOfCurves.Values);
                foreach (Dictionary<CharacterizedPrecursor, ElutionCurve> dic in lastDicOfCurves.Keys)
                    if (lastDicOfCurves[dic] < maxMed)
                        dicOfCurves.Add(dic, lastDicOfCurves[dic]);
                
                //Normalize already computed correlation factors for the remaning curves (sum must equal 1)
                double sumOfCorr = 0.0;
                foreach (Dictionary<CharacterizedPrecursor, ElutionCurve> dicOfCurve in dicOfCurves.Keys)
                    sumOfCorr += dicOfCorrelations[dicOfCurve];

                foreach (Dictionary<CharacterizedPrecursor, ElutionCurve> dicOfCurve in dicOfCurves.Keys)
                    dicOfCorrelations[dicOfCurve] /= sumOfCorr;

                //Compute average from weighted curves
                Dictionary<CharacterizedPrecursor, double> average = new Dictionary<CharacterizedPrecursor, double>();
                foreach (Dictionary<CharacterizedPrecursor, ElutionCurve> dicOfCurve in dicOfCurves.Keys)
                {
                    Dictionary<CharacterizedPrecursor, double> areas = GetAreas(dicOfCurve);
                    foreach (CharacterizedPrecursor cPep in areas.Keys)
                    {
                        if (!average.ContainsKey(cPep))
                            average.Add(cPep, 0);
                        average[cPep] += areas[cPep] * dicOfCorrelations[dicOfCurve];
                    }
                }

                //Compute correlation between average and curves
                List<double> corrs = new List<double>();
                foreach (Dictionary<CharacterizedPrecursor, ElutionCurve> dicOfCurve in dicOfCurves.Keys)
                {
                    Dictionary<CharacterizedPrecursor, double> elution = new Dictionary<CharacterizedPrecursor, double>();
                    foreach (CharacterizedPrecursor cPep in average.Keys)
                        if (dicOfCurve.ContainsKey(cPep))
                            elution.Add(cPep, dicOfCurve[cPep].Area);
                        else
                            elution.Add(cPep, 0);
                    double tmp =  1.0;
                    if(elution.Count > 1)
                        tmp = Math.Abs(MathNet.Numerics.Statistics.Correlation.Pearson(average.Values, elution.Values));

                    dicOfCorrelations[dicOfCurve] = tmp;
                    corrs.Add(tmp);
                }

                //Remove worst curves
                double medianCorr = MathNet.Numerics.Statistics.Statistics.Median(corrs);
                double maxCorr  = medianCorr + 0.5 * MathNet.Numerics.Statistics.Statistics.Variance(corrs);
                Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double> dicOfCurves2 = new Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double>();
                foreach (Dictionary<CharacterizedPrecursor, ElutionCurve> dic in dicOfCurves.Keys)
                    if (dicOfCorrelations[dic] < maxCorr)
                        dicOfCurves2.Add(dic, dicOfCurves[dic]);

                lastDicOfCurves = dicOfCurves2;
            }//End of While nbRun not exhausted

            Dictionary<CharacterizedPrecursor, double> bestAverage = new Dictionary<CharacterizedPrecursor, double>();
            foreach (Dictionary<CharacterizedPrecursor, ElutionCurve> dicOfCurve in lastDicOfCurves.Keys)
            {
                Dictionary<CharacterizedPrecursor, double> areas = GetAreas(dicOfCurve);
                foreach (CharacterizedPrecursor cPep in areas.Keys)
                {
                    if (!bestAverage.ContainsKey(cPep))
                        bestAverage.Add(cPep, 0);
                    bestAverage[cPep] += areas[cPep] * dicOfCorrelations[dicOfCurve];
                }
            }
            return bestAverage;
        }
        /*
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
        }//*/

        public static Dictionary<CharacterizedPrecursor, double> GetAreas(Dictionary<CharacterizedPrecursor, ElutionCurve> curves)
        {
            Dictionary<CharacterizedPrecursor, double> newCurves = new Dictionary<CharacterizedPrecursor, double>();
            foreach (CharacterizedPrecursor cPep in curves.Keys)
                newCurves.Add(cPep, curves[cPep].Area * cPep.PrecursorLossNormalizeFactor);// * cPep.FragmentLossNormalizeFactor);

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
                    if (charPrec.eCurve.Coefficients != null && charPrec.eCurve.Coefficients.Length == 3)
                        curveStr += "," + charPrec.eCurve.Coefficients[0] + "x^2 + " + charPrec.eCurve.Coefficients[1] + "x" + charPrec.eCurve.Coefficients[2];
                    else
                        curveStr += ",NA";
            writerCumul.AddLine(curveStr);

            foreach(double keyMz in characterizedPeptides.Keys)
            {
                //Get the list of precursors to characterize            
                foreach (Sample mixedSample in MixedSamples)
                {
                    Dictionary<double, MixedPrecursor> mixedPrecursors = GetMixedPrecursors(mixedSample, mixedResult, dbOptions, characterizedPeptides);

                    if (mixedPrecursors.ContainsKey(keyMz))
                    {
                        // Compute Max Flow for this precursor
                        Dictionary<CharacterizedPrecursor, double> ratios = GetRatiosFromMaxFlow(dbOptions, characterizedPeptides, mixedPrecursors[keyMz]);

                        string resultStr = vsCSV.GetFileName(mixedSample.sSDF) + "," + keyMz;
                        foreach (double precursor in characterizedPeptides.Keys)
                            foreach (CharacterizedPrecursor charPrec in characterizedPeptides[precursor].Values)
                                if (ratios.ContainsKey(charPrec))
                                    resultStr += "," + ratios[charPrec];
                                else
                                    resultStr += ",0";
                        writerCumul.AddLine(resultStr);

                        // Export results in a file
                        vsCSVWriter writerRatio = new vsCSVWriter(dbOptions.OutputFolder + vsCSV.GetFileName_NoExtension(mixedSample.sSDF) + "_" + keyMz + "MZ.csv");
                        writerRatio.AddLine("Sequence,Area");
                        writerRatio.AddLine("Total," + mixedPrecursors[keyMz].eCurve.Area);
                        foreach (CharacterizedPrecursor charPep in ratios.Keys)
                            writerRatio.AddLine(charPep.Peptide.Sequence + "," + ratios[charPep]);
                        writerRatio.WriteToFile();
                    }
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

                    //Don't try to characterize mixed precursors if there is less than three scans
                    if(mixedPrecursor.Queries.Count >= 3)
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
                            //Don't keep precursors if they are not well characterized (unfragmented or misasigned)
                            if (cPrec.Fragments.Count >= cPrec.Peptide.Length - 2)
                            {
                                if (!spikes.ContainsKey(mzKey))
                                    spikes.Add(mzKey, new Dictionary<Sample, CharacterizedPrecursor>());
                                if (!spikes[mzKey].ContainsKey(spikedSample))
                                    spikes[mzKey].Add(spikedSample, cPrec);
                                else
                                    Console.WriteLine("Twice??");
                            }
                        }
                    }
                }//End of foreach mzKey
            }//End of foreach spiked sample

            //Normalize intensities based on average area of each precursor
            foreach (double mzKey in spikes.Keys)
            {
                //foreach (CharacterizedPrecursor precursor in spikes[mzKey].Values)
                //    precursor.NormalizePrecursor(spikes[mzKey].Values);

                foreach(CharacterizedPrecursor precursor in spikes[mzKey].Values)
                    precursor.NormalizeFragments(spikes[mzKey].Values, dbOptions);
            }//*/
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
        
        public static Dictionary<CharacterizedPrecursor, double> GetRatiosFromMaxFlow(DBOptions dbOptions, Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> spikes, MixedPrecursor mixedPrecursor)        
        {
            int minProducts = 4;
            int maxProducts = 14;
            Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double> DicOfCurveErrors = new Dictionary<Dictionary<CharacterizedPrecursor, ElutionCurve>, double>();
            
            for (int nbProductsToKeep = minProducts; nbProductsToKeep <= maxProducts; nbProductsToKeep++)
            {
                int nbIgnoredSpectrum = 0;
                List<CharacterizedPrecursor> Isomers = new List<CharacterizedPrecursor>();
                foreach (double mz in spikes.Keys)
                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(mz, mixedPrecursor.MZ, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                        foreach (Sample sample in spikes[mz].Keys)
                            Isomers.Add(spikes[mz][sample]);

                Dictionary<double, int> dicOfCommonFragments = GetCommonFragmentMz(Isomers, nbProductsToKeep);
                Dictionary<CharacterizedPrecursor, List<ProductMatch>> matches = new Dictionary<CharacterizedPrecursor, List<ProductMatch>>();
                foreach (CharacterizedPrecursor cPrec in Isomers)
                    matches.Add(cPrec, cPrec.GetCombinedMatches(dicOfCommonFragments, dbOptions));

                //foreach(CharacterizedPrecursor precursor in Isomers)
                //    precursor.NormalizeFragments(Isomers, matches);                

                double cumulError = 0;
                Dictionary<CharacterizedPrecursor, ElutionCurve> curves = new Dictionary<CharacterizedPrecursor, ElutionCurve>();
                foreach (Query query in mixedPrecursor.Queries)
                {
                    double timeInMilliSeconds = query.spectrum.RetentionTimeInMin * 60.0 * 1000.0;
                    double overFlow = 0;
                    double underFlow = 0;
                    double percentError = 0;
                    Dictionary<CharacterizedPrecursor, MaxFlowResult> finalRatios = LaunchMaxFlowFromSpectrum(matches, 1000, query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                            0, query.spectrum.PrecursorIntensityPerMilliSecond * query.spectrum.InjectionTime, ref overFlow, ref underFlow, ref percentError, dbOptions.ConSole);

                    cumulError += percentError;
                    if (percentError < 0.5)
                    {
                        //UsedPrecursorArea += localArea;
                        //Dictionary<Sample, double> avgQuantifiedRatios = new Dictionary<Sample, double>();

                        foreach (CharacterizedPrecursor cPep in finalRatios.Keys)
                        {
                            if (!curves.ContainsKey(cPep))
                                curves.Add(cPep, new ElutionCurve());

                            curves[cPep].AddPoint(timeInMilliSeconds, finalRatios[cPep].Ratio * query.spectrum.PrecursorIntensityPerMilliSecond);
                        }
                    }
                    else
                        nbIgnoredSpectrum++;

                    if (nbIgnoredSpectrum * 2 > mixedPrecursor.Queries.Count)
                        break;
                }//End of foreach query
                
                //mixedPrecursor.PeptideRatios = curves;

                if (nbIgnoredSpectrum * 2 < mixedPrecursor.Queries.Count)
                {
                    if (nbIgnoredSpectrum > 0)
                        Console.WriteLine("Ignored Spectrum : " + nbIgnoredSpectrum + " / " + mixedPrecursor.Queries.Count);

                    foreach (CharacterizedPrecursor cPep in curves.Keys)
                        curves[cPep].Compute();

                    Dictionary<CharacterizedPrecursor, ElutionCurve> curvesToKeep = new Dictionary<CharacterizedPrecursor, ElutionCurve>();
                    foreach (CharacterizedPrecursor cPep in curves.Keys)
                        if (curves[cPep].Area > 0)
                            curvesToKeep.Add(cPep, curves[cPep]);

                    if (curvesToKeep.Count > 0)
                        DicOfCurveErrors.Add(curvesToKeep, cumulError);
                    /*
                        if (!values.ContainsKey(sRatio))
                            values.Add(sRatio, curves[sRatio].Area);
                        else
                            values[sRatio] += curves[sRatio].Area;
                    }//*/
                }
            }//End of for each nbProduct            

            Dictionary<CharacterizedPrecursor, double> averagedValues = mixedPrecursor.ComputePeptideRatios(DicOfCurveErrors);
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

        public static Dictionary<CharacterizedPrecursor, MaxFlowResult> LaunchMaxFlowFromSpectrum(Dictionary<CharacterizedPrecursor, List<ProductMatch>> ratiosToFit,  
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
            foreach (CharacterizedPrecursor key in ratiosToFit.Keys)
            {
                resultPerSample.Add(key, result[i]);
                i++;
            }
            return resultPerSample;
        }

        private static List<MaxFlowResult> GetResultList(List<double> solution, int precision, double underFlow, double sumOfIntensities)
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
    }
}
