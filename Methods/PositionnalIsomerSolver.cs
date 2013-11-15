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
    public class PositionnalIsomerSolver
    {
        public static Dictionary<Sample, double> GetNormalizationFactors(Result spikedResult, Result stableResult, DBOptions dbOptions, int charge,
                                int nbProductsMin, int nbProductsMax, bool smoothedPrecursor, int precision, int mflowReturnType, ref int nbProductsUsed)
        {
            Samples ProjectRatios = spikedResult.samples;
            Samples ProjectStable = stableResult.samples;
            Dictionary<Sample, double> tmpRatioNormalizer = new Dictionary<Sample, double>();
            foreach (Sample sample in ProjectRatios)
                tmpRatioNormalizer.Add(sample, 1.0);

            Dictionary<Sample, Dictionary<Sample, double>> tmpRatios = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectStable, tmpRatioNormalizer,
                                                            ref nbProductsUsed, mflowReturnType,
                                                            nbProductsMin, nbProductsMax, smoothedPrecursor,
                                                            precision, spikedResult, stableResult, charge);
            if (tmpRatios != null && tmpRatios.Count > 0)
            {
                Dictionary<Sample, double> RatioNormalizer = new Dictionary<Sample, double>();
                //foreach (Dictionary<double, double> listSample in tmpRatios.Values)
                foreach (Sample sampleMixed in tmpRatios.Keys)
                {
                    double sum = 0;
                    int nbUsedPeptides = 0;
                    foreach (double val in tmpRatios[sampleMixed].Values)
                    {
                        if (val > 0)
                            nbUsedPeptides++;

                        sum += val;
                    }
                    foreach (Sample sRatio in tmpRatios[sampleMixed].Keys)
                        if (!RatioNormalizer.ContainsKey(sRatio))
                            RatioNormalizer.Add(sRatio, 1.0);

                    double aim = 1.0 / (double)nbUsedPeptides;
                    foreach(Sample sRatio in tmpRatios[sampleMixed].Keys)
                        if (tmpRatios[sampleMixed][sRatio] > 0)
                        {
                            double aimRatio = aim / (tmpRatios[sampleMixed][sRatio] / sum);

                            //double aimRatioLog = Math.Log(aimRatio, 2);
                            //double aimRatioMult = Math.Pow(aimRatio, 2);

                            //Console.WriteLine("ADA " + aimRatioLog + "    " + aimRatioMult);
                            if (RatioNormalizer[sRatio] == 1.0 || Math.Abs(1.0 - aimRatio) < Math.Abs(1.0 - RatioNormalizer[sRatio]))
                                RatioNormalizer[sRatio] = aimRatio;
                        }
                }

                Dictionary<Sample, double> sampleNormalizer = new Dictionary<Sample, double>();
                foreach (Sample sampleMixed in tmpRatios.Keys)
                    foreach (Sample sampleRatio in tmpRatios[sampleMixed].Keys)
                    {
                        if(!sampleNormalizer.ContainsKey(sampleRatio))
                            sampleNormalizer.Add(sampleRatio, 0);

                        sampleNormalizer[sampleRatio] += tmpRatios[sampleMixed][sampleRatio] * RatioNormalizer[sampleRatio];
                    }
                /*
                double sumOfSamples = 0;
                foreach (Sample sample in sampleNormalizer.Keys)
                    sumOfSamples += sampleNormalizer[sample];
                double avg = sumOfSamples / (double) sampleNormalizer.Count;

                foreach (Sample sample in sampleNormalizer.Keys)
                    RatioNormalizer[sample] /= avg / sampleNormalizer[sample];
                //*/
                return RatioNormalizer;
            }
            else
            {
                return null;
            }
        }

        public static void Solve(Samples ProjectRatios, Samples ProjectStable, Samples ProjectMixed, DBOptions dbOptions,
                                int nbProductsMin = 4, int nbProductsMax = 7, bool smoothedPrecursor = false, int precision = 1000, int maxCharge = 4)
        {
            int mflowReturnType = 0;

            //Precompute Spiked peptide identifications
            Result spikedResult = Propheus.Start(dbOptions, ProjectRatios, false, false, false);

            Result stableResult = null;
            if(ProjectStable != null && ProjectStable.Count > 0)
                stableResult = Propheus.Start(dbOptions, ProjectStable, false, false, false);

            Result mixedResult = Propheus.Start(dbOptions, ProjectMixed, false, false, false);
                        
            for (int charge = 2; charge <= maxCharge; charge++)
            {
                //Get normalization factor from the stable mix
                int nbProductsUsed = 0;

                Dictionary<Sample, double> RatioNormalizer = null;
                if (ProjectStable != null && ProjectStable.Count > 0)
                {
                    RatioNormalizer = GetNormalizationFactors(spikedResult, stableResult, dbOptions, charge, nbProductsMin, nbProductsMax,
                                                              smoothedPrecursor, precision, mflowReturnType, ref nbProductsUsed);
                    if (nbProductsUsed > 0)
                    {
                        nbProductsMin = nbProductsUsed;
                        nbProductsMax = nbProductsUsed;
                    }
                }
                if(RatioNormalizer == null)
                {
                    RatioNormalizer = new Dictionary<Sample, double>();
                    foreach(Sample sRatio in ProjectRatios)
                        RatioNormalizer.Add(sRatio, 1.0);
                }

                //Solve positional isomers in the mixed spectrums
                Dictionary<Sample, Dictionary<Sample, double>> ratios = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectMixed, RatioNormalizer, ref nbProductsUsed, mflowReturnType,
                                                            nbProductsMin, nbProductsMax, smoothedPrecursor, precision, spikedResult, mixedResult, charge);

                if (ratios == null || ratios.Count == 0)
                    dbOptions.ConSole.WriteLine("No precursor found with charge " + charge);
                else
                {
                    //Export results to a csv file
                    vsCSVWriter writerCumul = new vsCSVWriter(dbOptions.OutputFolder + "CumulRatios9_charge" + charge + ".csv");
                    foreach(Sample mixedSample in ProjectMixed)
                    {
                        string lineCumulRatio = mixedSample.sSDF;
                        foreach (Sample ratioSample in ProjectRatios)
                            if (ratios[mixedSample].ContainsKey(ratioSample))
                                lineCumulRatio += "," + ratios[mixedSample][ratioSample];
                            else
                                lineCumulRatio += ",0.0";

                        writerCumul.AddLine(lineCumulRatio);
                    }
                    writerCumul.WriteToFile();

                    //Export Statistics of the approach
                    vsCSVWriter writerStats = new vsCSVWriter(dbOptions.OutputFolder + "Stats_Charge" + charge + ".txt");
                    writerStats.AddLine("ProjectRatios : " + ProjectRatios.FileName);
                    writerStats.AddLine("ProjectStable : " + (ProjectStable == null ? "none" : ProjectStable.FileName));
                    writerStats.AddLine("ProjectMixed : " + ProjectMixed.FileName);
                    writerStats.AddLine("Nb Of products considered : " + nbProductsUsed);
                    writerStats.AddLine("Precision : " + precision);
                    writerStats.AddLine("mflowReturnType : " + mflowReturnType);
                    writerStats.AddLine("smoothedPrecursor : " + smoothedPrecursor);
                    //writerStats.AddLine("aim (for stable ratios) : " + aim);
                    foreach(Sample sample in RatioNormalizer.Keys)
                        writerStats.AddLine("Normalization for " + vsCSV.GetFileName_NoExtension(sample.sSDF) + " : " + RatioNormalizer[sample]);

                    //Compute possible contamination of original spiked peptides
                    writerStats.AddLine(" -- Cross Talk between Spiked Samples -- ");
                    Dictionary<Sample, Dictionary<Sample, double>> ratiosForCrossTalk = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectRatios, RatioNormalizer, ref nbProductsUsed, mflowReturnType,
                                                                            nbProductsMin, nbProductsMax, smoothedPrecursor, precision, spikedResult, spikedResult, charge);

                    //Export to the same statistic csv file
                    foreach(Sample sample in ratiosForCrossTalk.Keys)
                    {
                        writerStats.AddLine(vsCSV.GetFileName_NoExtension(sample.sSDF));
                        double sum = 0;
                        foreach (double val in ratiosForCrossTalk[sample].Values)
                            sum += val;
                        writerStats.AddLine("Coverage : " + (float)(100.0 * ratiosForCrossTalk[sample][sample] / sum));
                        string line = "Spreading : ";
                        foreach (double val in ratiosForCrossTalk[sample].Values)
                            line += (val / sum).ToString() + ",";
                        writerStats.AddLine(line);
                    }
                    writerStats.WriteToFile();
                }
            }
            //Export results of searches
            spikedResult.Export(0.05);
            mixedResult.Export(0.05);
        }    
    
        public static Dictionary<Sample, Dictionary<Sample, double>> ComputeMaxFlows(
                                DBOptions dbOptions,
                                Samples ProjectRatios,
                                Samples ProjectMixed,
                                Dictionary<Sample, double> RatioNormalizer,
                                ref int nbProductsUsed,
                                int mflowReturnType,
                                int nbProductMin,// 5,
                                int nbProductMax,// 5,
                                bool smoothedPrecursor,
                                int precision,
                                Result precomputedSpiked,
                                Result precomputedMixed,
                                int chargeToConsider,
                                double aim4StableRatio = -1)
        {
            Dictionary<Sample, double> peptideMasses = new Dictionary<Sample, double>();
            foreach (Sample sample in ProjectRatios)
                peptideMasses.Add(sample, Peptide.ComputeMonoisotopicMass(sample.nameColumn));//*/

            //Get PSMs
            Result mixedResult = precomputedMixed;
            Result spikedResult = precomputedSpiked;

            Dictionary<int, Dictionary<double, Dictionary<Sample, List<ProductMatch>>>> DicOfRatios = new Dictionary<int, Dictionary<double, Dictionary<Sample, List<ProductMatch>>>>();
            //Dictionary<int, Dictionary<double, List<double>>>               DicOfPrecursorAreas = new Dictionary<int, Dictionary<double, List<double>>>();
            Dictionary<int, Dictionary<double, Dictionary<Sample, double>>> DicOfNormalizeFactor = new Dictionary<int, Dictionary<double, Dictionary<Sample, double>>>();
            Dictionary<int, double>                                         DicOfErrors = new Dictionary<int,double>();
            Dictionary<int, List<double>>                                   DicOfFragmentMz = new Dictionary<int, List<double>>();

            for (int nbProductsToKeep = nbProductMin; nbProductsToKeep <= nbProductMax; nbProductsToKeep++)
            {
                Dictionary<double, Dictionary<Sample, List<ProductMatch>>> ratios = new Dictionary<double, Dictionary<Sample, List<ProductMatch>>>();
                Dictionary<double, Dictionary<Sample, double>> PrecursorAreas = new Dictionary<double, Dictionary<Sample, double>>();
                Dictionary<double, Dictionary<Sample, double>> NormalizeFactor = new Dictionary<double, Dictionary<Sample, double>>();
                List<double> ListFragmentMz = new List<double>();

                BuildSinglePeptideVirtualSpectrum(spikedResult, smoothedPrecursor, nbProductsToKeep, RatioNormalizer,
                                                        ref ratios, ref PrecursorAreas, ref NormalizeFactor, ref ListFragmentMz, chargeToConsider);

                DicOfRatios.Add(nbProductsToKeep, ratios);
                //DicOfPrecursorAreas.Add(nbProductsToKeep, PrecursorAreas);
                DicOfNormalizeFactor.Add(nbProductsToKeep, NormalizeFactor);
                DicOfFragmentMz.Add(nbProductsToKeep, ListFragmentMz);
                DicOfErrors.Add(nbProductsToKeep, 0);
            }

            double smallestError = double.MaxValue;
            Dictionary<Sample, Dictionary<Sample, double>> bestListOfSumOfRatio = null;
            Dictionary<Sample, List<Dictionary<int, string>>> bestDicOfResults = null;
            for (int nbProductsToKeep = nbProductMin; nbProductsToKeep <= nbProductMax; nbProductsToKeep++)
            {
                long iterError = 0;
                double cumulPercentError = 0;
                Dictionary<Sample, Dictionary<Sample, double>> listOfSumOfRatio = new Dictionary<Sample, Dictionary<Sample, double>>();
                Dictionary<Sample, List<Dictionary<int, string>>> dicOfResultsPerSample = new Dictionary<Sample, List<Dictionary<int, string>>>();
                foreach (Sample mixedSample in ProjectMixed)
                {
                    //Dictionary<double, double> listOfSumOfRatioPerPrecursor = new Dictionary<double, double>();
                    Dictionary<Sample, double> DicOfSumOfRatioPerPeptide = new Dictionary<Sample, double>();
                    dicOfResultsPerSample.Add(mixedSample, new List<Dictionary<int, string>>());
                    //Dictionary<ProductSpectrum, bool> doneSpectrum = new Dictionary<ProductSpectrum, bool>();
                    //Dictionary<PeptideSpectrumMatch, ProductSpectrum> DicOfSpectrum = new Dictionary<ProductSpectrum, PeptideSpectrumMatch>();
                    Dictionary<ProductSpectrum, PeptideSpectrumMatch> DicOfSpectrum = new Dictionary<ProductSpectrum, PeptideSpectrumMatch>();

                    //Get groups of spectrum (split spectrum of different masses) into different groups
                    Dictionary<double, List<ProductSpectrum>> DicOfSpectrumMasses = new Dictionary<double, List<ProductSpectrum>>();
                    foreach (Query query in mixedResult.queries)
                    {
                        if (query.sample == mixedSample && !DicOfSpectrum.ContainsKey(query.spectrum) && query.precursor.Charge == chargeToConsider)
                        {                            
                            PeptideSpectrumMatch psmToKeep = null;
                            foreach (PeptideSpectrumMatch psm in query.psms)
                                foreach(Sample sRatio in ProjectRatios)
                                    if (sRatio.nameColumn.CompareTo(psm.Peptide.Sequence) == 0)
                                        psmToKeep = psm;
                            
                            if (psmToKeep != null)
                            {
                                double foundKey = 0.0;
                                foreach (double key in DicOfRatios[nbProductsToKeep].Keys)
                                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(psmToKeep.Peptide.MonoisotopicMass, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                        foundKey = key;
                                if(foundKey == 0.0)
                                    dbOptions.ConSole.WriteLine("Mass key not found.... not good!");

                                if(!DicOfSpectrumMasses.ContainsKey(foundKey))
                                {
                                    List<ProductSpectrum> listOfSpectrum = new List<ProductSpectrum>();
                                    listOfSpectrum.Add(query.spectrum);
                                    DicOfSpectrumMasses.Add(foundKey, listOfSpectrum);
                                }
                                else
                                {
                                    DicOfSpectrumMasses[foundKey].Add(query.spectrum);
                                }
                                DicOfSpectrum.Add(query.spectrum, psmToKeep);
                            }
                        }
                    }
                    
                    //Cycle through groups of spectrum and evaluate metrics as if they where different samples
                    foreach (double key in DicOfSpectrumMasses.Keys)
                    {
                        Dictionary<Sample, double> sumOfRatio = new Dictionary<Sample, double>();

                        List<ProductSpectrum> sortedSpectrum = DicOfSpectrumMasses[key];
                        sortedSpectrum.Sort(ProductSpectrum.AscendingRetentionTimeComparison);
                        PeptideSpectrumMatches psmMatches = new PeptideSpectrumMatches();
                        foreach (ProductSpectrum ps in sortedSpectrum)
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(DicOfSpectrum[ps].Peptide.MonoisotopicMass, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                psmMatches.Add(DicOfSpectrum[ps]);

                        psmMatches.Sort(PeptideSpectrumMatches.AscendingRetentionTime);

                        // -- Export best spectrum -- //
                        PeptideSpectrumMatch bestPsm = null;
                        foreach (PeptideSpectrumMatch psm in psmMatches)
                            if (bestPsm == null || psm.Query.spectrum.TotalIntensity > bestPsm.Query.spectrum.TotalIntensity)
                                bestPsm = psm;
                        mixedResult.ExportFragments(bestPsm);
                        // -- ******************** -- //


                        //Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor = psmMatches.ComputeMsMsNormalizationFactors();
                        
                        Dictionary<int, string> dicOfResults = new Dictionary<int, string>();
                        
                        double UsedPrecursorArea = 0.0;                        
                        double LastTimeStamp = 0;
                        //double TotalElapsedTime = 0;
                        //foreach(ProductSpectrum key in sortedSpectrum)
                        foreach (PeptideSpectrumMatch psm in psmMatches)
                        {
                            double overFlow = 0;
                            double underFlow = 0;
                            double percentError = 0;
                            Dictionary<Sample, List<ProductMatch>> matches = new Dictionary<Sample, List<ProductMatch>>();

                            //for (int i = 0; i < peptideMasses.Count; i++)
                            foreach(Sample sRatio in peptideMasses.Keys)
                                if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(psm.Peptide.MonoisotopicMass, peptideMasses[sRatio], dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                    matches.Add(sRatio, DicOfRatios[nbProductsToKeep][key][sRatio]);
                            Dictionary<Sample, double> finalRatios = LaunchMaxFlowFromSpectrum(matches, precision, psm.Query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                                    mflowReturnType, psm.Query.spectrum.PrecursorIntensityPerMilliSecond * psm.Query.spectrum.InjectionTime, ref overFlow, ref underFlow, ref percentError, dbOptions.ConSole);

                            //double normSumRatio = 0;
                            //for (int i = 0; i < finalRatios.Count; i++)
                            //    normSumRatio += finalRatios[i] * DicOfNormalizeFactor[nbProductsToKeep][key][i];

                            //double sumRatio = 0.0;
                            //for (int i = 0; i < finalRatios.Count; i++)
                            //    sumRatio += finalRatios[i];

                            double ElapsedTime = (psm.Query.spectrum.RetentionTimeInMin - LastTimeStamp) * 60.0 * 1000.0;

                            double localArea = psm.Query.spectrum.PrecursorIntensityPerMilliSecond * ElapsedTime;                            

                            if (LastTimeStamp > 0)
                            {
                                //if (IntensityPerUnitOfTime > 1.75 * LastIntensityPerUnitOrTime)
                                //    Console.WriteLine("oops?");
                                //precursorArea += localArea;
                                if (percentError < 0.5)//25)// && !double.IsNaN(normSumRatio) && normSumRatio > 0)//0.95
                                {
                                    UsedPrecursorArea += localArea;
                                    Dictionary<Sample, double> avgQuantifiedRatios = new Dictionary<Sample, double>();

                                    foreach (Sample sRatio in finalRatios.Keys)
                                    {
                                        //double avgRatio = (finalRatios[i] / sumRatio) * DicOfNormalizeFactor[nbProductsToKeep][key][i] * localArea;
                                        double avgRatio = finalRatios[sRatio] * DicOfNormalizeFactor[nbProductsToKeep][key][sRatio] * localArea * RatioNormalizer[sRatio];
                                        if (double.IsNaN(avgRatio))
                                            dbOptions.ConSole.WriteLine("Oops, NaN in ratios");
                                        avgQuantifiedRatios.Add(sRatio, avgRatio);
                                    }

                                    string strRatios = psm.Query.spectrum.RetentionTimeInMin.ToString() + "," + psm.Query.spectrum.PrecursorIntensityPerMilliSecond;

                                    foreach (Sample sRatio in ProjectRatios)
                                    {
                                        if (avgQuantifiedRatios.ContainsKey(sRatio))
                                        {
                                            if (!sumOfRatio.ContainsKey(sRatio))
                                                sumOfRatio.Add(sRatio, 0.0);
                                            sumOfRatio[sRatio] += avgQuantifiedRatios[sRatio];
                                            strRatios += "," + finalRatios[sRatio] * DicOfNormalizeFactor[nbProductsToKeep][key][sRatio] * RatioNormalizer[sRatio] * psm.Query.spectrum.PrecursorIntensityPerMilliSecond;
                                        }
                                        else
                                            strRatios += ",0.0";
                                    }
                                    foreach (double mz in DicOfFragmentMz[nbProductsToKeep])
                                    {
                                        double intCumul = 0;
                                        foreach (MsMsPeak peak in psm.Query.spectrum.Peaks)
                                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(mz, peak.MZ, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                                intCumul += peak.Intensity;
                                        strRatios += "," + intCumul;
                                    }

                                    dicOfResults.Add(psm.Query.spectrum.ScanNumber, strRatios);
                                    //TotalElapsedTime += ElapsedTime;                                                                                
                                }
                                else
                                    dbOptions.ConSole.WriteLine("Bad MaxFlow computation : " + percentError);

                                //foreach (double ratio in finalRatios.Values)
                                //    if (ratio == 0)
                                //        percentError += 1.0 / (double)finalRatios.Count;//*/
                                cumulPercentError += percentError;
                                iterError++;
                            }

                            LastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;    
                        }//end of foreach psm

                        //Interpolate unused spectrum
                        double ComputedPrecursorArea = psmMatches.ComputePrecursorArea(smoothedPrecursor);
                        double factorOfUnusedLocalAreas = (ComputedPrecursorArea - UsedPrecursorArea) / UsedPrecursorArea;

                        dicOfResultsPerSample[mixedSample].Add(dicOfResults);

                        //Use SumOfRatio as the fraction for the entire Precursor area
                        /*double sumArea = 0;
                        foreach (double val in sumOfRatio)
                            sumArea += val;
                        
                        List<double> relativeSumOfRatio = new List<double>();
                        for (int i = 0; i < sumOfRatio.Count; i++)
                            sumOfRatio[i] = (sumOfRatio[i] / sumArea) * precursorArea;
                        //*/
                        
                        foreach (Sample ratioSample in sumOfRatio.Keys)
                        {
                            if(!DicOfSumOfRatioPerPeptide.ContainsKey(ratioSample))
                                DicOfSumOfRatioPerPeptide.Add(ratioSample, 0);
                            DicOfSumOfRatioPerPeptide[ratioSample] += sumOfRatio[ratioSample] + sumOfRatio[ratioSample] * factorOfUnusedLocalAreas;
                        }
                    }//end of foreach precursor mass

                    listOfSumOfRatio.Add(mixedSample, DicOfSumOfRatioPerPeptide);//Per Precursor

                }//end of foreach mixed sample
                /*
                cumulPercentError /= (double)iterError;
                if (aim4StableRatio > 0)
                {
                    foreach (Dictionary<Sample,double> dicOfSum in listOfSumOfRatio.Values)
                    {
                        double sum = 0;
                        foreach (double val in dicOfSum.Values)
                            sum += val;

                        foreach (double val in dicOfSum.Values)
                            cumulPercentError += Math.Abs(aim4StableRatio - val / sum);
                    }
                }//*/
                
                //cumulPercentError /= (double)nbProductsToKeep;//Average error, per peak/fragment
                if (cumulPercentError < smallestError)
                {
                    nbProductsUsed = nbProductsToKeep;
                    smallestError = cumulPercentError;
                    bestListOfSumOfRatio = listOfSumOfRatio;
                    bestDicOfResults = dicOfResultsPerSample;
                }
            }//end of foreach nbProductUsed

            if (bestDicOfResults != null)
            {
                dbOptions.ConSole.WriteLine("Best Number of products : " + nbProductsUsed);
                foreach (Sample mixedSample in bestDicOfResults.Keys)
                {//Dictionary<Sample, List<Dictionary<int, string>>> 
                    vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + mixedSample.nameColumn + "_Ratios7_Charge" + chargeToConsider + ".csv");
                    string title = "Retention Time,Precursor Intensity";
                    for(int i = 0; i < ProjectRatios.Count; i++)
                        title += "," + ProjectRatios[i].nameColumn;
                    
                    foreach (double mz in DicOfFragmentMz[nbProductsUsed])
                        title += "," + mz;
                    writer.AddLine(title);
                    foreach (Dictionary<int, string> dic in bestDicOfResults[mixedSample])
                    {
                        List<int> scans = new List<int>(dic.Keys);
                        scans.Sort();
                        foreach (int scan in scans)
                            writer.AddLine(dic[scan]);
                    }
                    writer.WriteToFile();
                }
                //foreach (List<double> sampleRatio in bestListOfSumOfRatio)
                //    for (int i = 0; i < sampleRatio.Count; i++)
                //        sampleRatio[i] *= RatioNormalizer[i];
            }
            return bestListOfSumOfRatio;
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
            int bestIndex = 0;

            int iterSize = 1;
            double bestOverallError = double.MaxValue;
            List<long> bestLocalFlows = new List<long>();

            while (overError > 1 && iterSize < 10000)//anything less than 1 is an acceptable solution
            {
                bestIndex = -1;
                //double smallestUnderError = double.MaxValue;
                //double smallestOverError = overError;
                double worstFlowRate = 0.0;
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
                        if (double.IsNaN(tmpFlowRate) || double.IsInfinity(tmpFlowRate))
                            ConSole.WriteLine("schnit");
                        if(tmpFlowRate > worstFlowRate)
                        //if (tmpErrorPlus < overError && (tmpErrorMinus < smallestUnderError
                        //    || (tmpErrorMinus == smallestUnderError && tmpErrorPlus < smallestOverError)))
                        {
                            worstFlowRate = tmpFlowRate;
                            //smallestOverError = tmpErrorPlus;
                            //smallestUnderError = tmpErrorMinus;
                            bestIndex = i;
                        }
                        localFlows[i] += iterSize;
                    }
                }
                if (bestIndex != -1)
                {
                    localFlows[bestIndex] -= iterSize;
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
            }
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
        }//*/

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
                double intensity = peak.Intensity;// +peak.Intensity * boostRatioForSpectrum;
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
                    //newPm.obsIntensity *= sumOfProducts;
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
                    result = GetResult(solutions[0], precision, underFlow, sumOfIntensities);
                    break;
                case 1:
                    result = GetResult(solutions[1], precision, underFlow, sumOfIntensities);
                    break;
                case 2:
                    List<double> tmpAverage = new List<double>();                    
                    foreach (double val in average)
                        tmpAverage.Add(val);
                    result = GetResult(tmpAverage, precision, underFlow, sumOfIntensities);
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

        private static List<double> GetResult(List<double> solution, int precision, double underFlow, double sumOfIntensities)
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

        private static void BuildSinglePeptideVirtualSpectrum(Result precomputedResults, bool smoothPrecursor, 
                    int nbProductsToKeep, Dictionary<Sample, double> RatioNormalizer, ref Dictionary<double, Dictionary<Sample, List<ProductMatch>>> FinalSpikedProducts, 
                    ref Dictionary<double, Dictionary<Sample, double>> PrecursorAreas, ref Dictionary<double, Dictionary<Sample, double>> Normalizor, ref List<double> FragmentMz, int charge)
        {
            DBOptions dbOptions = precomputedResults.dbOptions;
            Samples Project = precomputedResults.samples;

            foreach (Sample sample in Project)
            {
                double peptideMass = Peptide.ComputeMonoisotopicMass(sample.nameColumn);

                double foundKey = 0.0;
                foreach (double key in FinalSpikedProducts.Keys)
                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(peptideMass, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                        foundKey = key;

                if (foundKey == 0.0)
                {
                    foundKey = peptideMass;
                    FinalSpikedProducts.Add(foundKey, new Dictionary<Sample,List<ProductMatch>>());
                    PrecursorAreas.Add(foundKey, new Dictionary<Sample, double>());
                    Normalizor.Add(foundKey, new Dictionary<Sample, double>());
                }
            }
            
            Dictionary<double, int> AllFragments = new Dictionary<double,int>();
            vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + "FragmentsUsed_" + nbProductsToKeep + "Products.csv");
            writer.AddLine("Precursor Mass, Peptide Sequence,Fragment,Pos,Charge,Mz,Intensity,Normalized Intensity");
            foreach (double foundKey in FinalSpikedProducts.Keys)
            {
                Peptide peptide = null;
                Dictionary<Sample, List<ProductMatch>> SpikedProducts = new Dictionary<Sample,List<ProductMatch>>();
                //List<double> PeakAvgIntensities = new List<double>();
                //Get list of PSMs per sample
                double averagePrecursorArea = 0;
                //double avgMsMsProductIntensity = 0;
                //int nbAveragedSpectrum = 0;
                int nbAverageArea = 0;
                Dictionary<Sample, PeptideSpectrumMatches> listOfPSMs = new Dictionary<Sample, PeptideSpectrumMatches>();
                foreach (Sample sample in Project)
                {
                    //double avgPeakIntensity = 0;
                    PeptideSpectrumMatches psmList = new PeptideSpectrumMatches();
                    foreach (Query query in precomputedResults.queries)
                    {
                        if (query.precursor.Charge == charge && query.precursor.sample == sample)
                        {
                            foreach (PeptideSpectrumMatch psm in query.psms)
                            {
                                if (sample.nameColumn.CompareTo(psm.Peptide.Sequence) == 0)// && psm.MatchingProducts > 4)//.ProbabilityScore() > 0.020)// && psm.ProbabilityScore() > bestPsmScore)
                                {
                                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(psm.Peptide.MonoisotopicMass, foundKey, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                    {
                                        peptide = psm.Peptide;
                                        psmList.Add(psm);
                                        //avgPeakIntensity += query.spectrum.PrecursorIntensity;
                                        //avgMsMsProductIntensity += query.spectrum.TotalIntensity;
                                        //nbAveragedSpectrum++;
                                    }
                                }
                            }
                        }
                    }

                    if (psmList.Count > 0)
                    {
                        //avgPeakIntensity /= (double)psmList.Count;
                        psmList.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
                        listOfPSMs.Add(sample, psmList);
                        List<ProductMatch> productList = psmList.GetCombinedSpectrum(precomputedResults.dbOptions, peptide, charge);

                        productList.Sort(ProductMatch.AscendingWeightComparison);

                        double precursorArea = psmList.ComputePrecursorArea(smoothPrecursor);
                        averagePrecursorArea += precursorArea;
                        nbAverageArea++;

                        SpikedProducts.Add(sample, productList);
                        PrecursorAreas[foundKey].Add(sample, precursorArea);
                        //PeakAvgIntensities.Add(avgPeakIntensity);
                    }
                    else
                    {
                        listOfPSMs.Add(sample, psmList);
                        PrecursorAreas[foundKey].Add(sample, 0);
                        //PeakAvgIntensities.Add(0);
                    }
                }
                //avgMsMsProductIntensity /= nbAveragedSpectrum;
                averagePrecursorArea /= (double)nbAverageArea;

                Dictionary<double, int> DicOfFragmentsToKeep = new Dictionary<double, int>();
                //Get List of desired fragments, and keep masses
                foreach(Sample sample in SpikedProducts.Keys)
                {
                    for (int j = SpikedProducts[sample].Count - 1; j > 0 && j >= SpikedProducts[sample].Count - nbProductsToKeep; j--)
                    {
                        if (!DicOfFragmentsToKeep.ContainsKey(SpikedProducts[sample][j].theoMz))
                            DicOfFragmentsToKeep.Add(SpikedProducts[sample][j].theoMz, 0);

                        DicOfFragmentsToKeep[SpikedProducts[sample][j].theoMz]++;
                    }
                }

                foreach (Sample sample in SpikedProducts.Keys)
                {
                    List<ProductMatch> list = listOfPSMs[sample].GetCombinedSpectrum(precomputedResults.dbOptions, peptide, charge, DicOfFragmentsToKeep);

                    foreach (ProductMatch match in list)
                        writer.AddLine(foundKey + "," + sample.nameColumn + "," + match.fragment + "," + match.fragmentPos + "," + match.charge + "," + match.theoMz + "," + match.obsIntensity + "," + match.normalizedIntensity);

                    FinalSpikedProducts[foundKey].Add(sample, list);
                }

                foreach(double key in DicOfFragmentsToKeep.Keys)
                    if(!AllFragments.ContainsKey(key))
                        AllFragments.Add(key, 1);
                    else
                        AllFragments[key]++;
                /*
                double avgPrecursors = 0;
                int nbInt = 0;
                foreach (double avg in PeakAvgIntensities)
                {
                    if (avg > 0)
                        nbInt++;
                    avgPrecursors += avg;
                }
                avgPrecursors /= (double)nbInt;
                /*
                for (int i = 0; i < Project.Count; i++)
                {
                    //Normalize each spectrum based on average precursor intensity
                    if (PeakAvgIntensities[i] > 0)
                    {
                        double y = (avgPrecursors - PeakAvgIntensities[i]) / PeakAvgIntensities[i];
                        for (int j = 0; j < FinalSpikedProducts[foundKey][i].Count; j++)
                            FinalSpikedProducts[foundKey][i][j].obsIntensity += y * FinalSpikedProducts[foundKey][i][j].obsIntensity;
                    }
                }//*/

                Dictionary<Sample, double> listOfsumOfProducts = new Dictionary<Sample, double>();
                double avgRatioSum = 0;
                int nbInt = 0;
                foreach (Sample sample in FinalSpikedProducts[foundKey].Keys)
                {
                    double sumOfProducts = 0;
                    for (int j = 0; j < FinalSpikedProducts[foundKey][sample].Count; j++)
                        sumOfProducts += FinalSpikedProducts[foundKey][sample][j].normalizedIntensity;
                    listOfsumOfProducts.Add(sample, sumOfProducts);

                    avgRatioSum += sumOfProducts; 

                    if (sumOfProducts > 0)
                        nbInt++;
                }
                avgRatioSum /= (double)nbInt;

                foreach (Sample sample in listOfsumOfProducts.Keys)
                {
                    if (listOfsumOfProducts[sample] > 0)
                    {
                        double lossOfPrecIntensity = Math.Log(averagePrecursorArea, 2) / Math.Log(PrecursorAreas[foundKey][sample], 2);
                        double lossofFragmentIntensity = Math.Pow(avgRatioSum, 2) / Math.Pow(listOfsumOfProducts[sample], 2);

                        //Normalizor[foundKey].Add((spectrumWeight / listOfsumOfProducts[i]) * (averagePrecursorArea / PrecursorAreas[foundKey][i]) * RatioNormalizer[i]);
                        //Normalizor[foundKey].Add(lossOfPrecIntensity * lossofFragmentIntensity);//(averagePrecursorArea / PrecursorAreas[foundKey][i]) * RatioNormalizer[i]);
                        Normalizor[foundKey].Add(sample, lossofFragmentIntensity * lossOfPrecIntensity);//(averagePrecursorArea / PrecursorAreas[foundKey][i]) * RatioNormalizer[i]);
                        //for (int j = 0; j < FinalSpikedProducts[foundKey][i].Count; j++)
                        //    FinalSpikedProducts[foundKey][i][j].normalizedIntensity /= listOfsumOfProducts[i];
                    }
                    else
                        Normalizor[foundKey].Add(sample, 1.0);
                }
            }

            foreach(double key in AllFragments.Keys)
                FragmentMz.Add(key);

            writer.WriteToFile();
        }
    }
}
