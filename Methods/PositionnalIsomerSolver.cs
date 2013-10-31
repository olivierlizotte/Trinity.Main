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
        public static List<double> GetNormalizationFactors(Result spikedResult, Result stableResult, DBOptions dbOptions, int charge,
                                int nbProductsMin, int nbProductsMax, bool smoothedPrecursor, int precision, int mflowReturnType, double aim, ref int nbProductsUsed)
        {
            Samples ProjectRatios = spikedResult.samples;
            Samples ProjectStable = stableResult.samples;
            List<double> tmpRatioNormalizer = new List<double>();
            foreach (Sample sample in ProjectRatios)
                tmpRatioNormalizer.Add(1.0);

            List<List<double>> tmpRatios = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectStable, tmpRatioNormalizer,
                                                            ref nbProductsUsed, mflowReturnType,
                                                            nbProductsMin, nbProductsMax, smoothedPrecursor,
                                                            precision, spikedResult, stableResult, charge, aim);
            if (tmpRatios != null && tmpRatios.Count > 0)
            {
                double sum = 0;
                foreach (double val in tmpRatios[0])
                    sum += val;

                List<double> RatioNormalizer = new List<double>();
                foreach (double val in tmpRatios[0])
                    RatioNormalizer.Add(aim / (val / sum));
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
            double aim = 1.0 / (double)ProjectRatios.Count;

            //Precompute Spiked peptide identifications
            Result spikedResult = Propheus.Start(dbOptions, ProjectRatios, false, false, false);

            Result stableResult = null;
            if(ProjectStable != null && ProjectStable.Count > 0)
                stableResult = Propheus.Start(dbOptions, ProjectStable, false, false, false);

            Result mixedResult = Propheus.Start(dbOptions, ProjectMixed, false, false, false);

            for (int charge = 2; charge < maxCharge; charge++)
            {
                //Get normalization factor from the stable mix
                int nbProductsUsed = 0;

                List<double> RatioNormalizer = null;
                if (ProjectStable != null && ProjectStable.Count > 0)
                {
                    RatioNormalizer = GetNormalizationFactors(spikedResult, stableResult, dbOptions, charge, nbProductsMin, nbProductsMax,
                                                              smoothedPrecursor, precision, mflowReturnType, aim, ref nbProductsUsed);
                    if (nbProductsUsed > 0)
                    {
                        nbProductsMin = nbProductsUsed;
                        nbProductsMax = nbProductsUsed;
                    }
                }
                if(RatioNormalizer == null)
                {
                    RatioNormalizer = new List<double>();
                    for (int i = 0; i < ProjectRatios.Count; i++)
                        RatioNormalizer.Add(1.0);
                }

                //Solve positional isomers in the mixed spectrums
                List<List<double>> ratios = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectMixed, RatioNormalizer, ref nbProductsUsed, mflowReturnType,
                                                            nbProductsMin, nbProductsMax, smoothedPrecursor, precision, spikedResult, mixedResult, charge);

                if (ratios == null || ratios.Count == 0)
                    Console.WriteLine("No precursor found with charge " + charge);
                else
                {
                    //Export results to a csv file
                    vsCSVWriter writerCumul = new vsCSVWriter(dbOptions.OutputFolder + "CumulRatios8_charge" + charge + ".csv");
                    for (int i = 0; i < ProjectMixed.Count; i++)
                    {
                        string lineCumulRatio = ProjectMixed[i].sSDF;
                        for (int k = 0; k < ratios[i].Count; k++)
                            lineCumulRatio += "," + ratios[i][k];
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
                    writerStats.AddLine("aim (for stable ratios) : " + aim);
                    for (int i = 0; i < RatioNormalizer.Count; i++)
                        writerStats.AddLine("Normalization for " + vsCSV.GetFileName_NoExtension(ProjectRatios[i].sSDF) + " : " + RatioNormalizer[i]);

                    //Compute possible contamination of original spiked peptides
                    writerStats.AddLine(" -- Cross Talk between Spiked Samples -- ");
                    List<List<double>> ratiosForCrossTalk = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectRatios, RatioNormalizer, ref nbProductsUsed, mflowReturnType,
                                                                            nbProductsMin, nbProductsMax, smoothedPrecursor, precision, spikedResult, spikedResult, charge);

                    //Export to the same statistic csv file
                    for (int i = 0; i < ratiosForCrossTalk.Count; i++)
                    {
                        writerStats.AddLine(vsCSV.GetFileName_NoExtension(ProjectRatios[i].sSDF));
                        double sum = 0;
                        foreach (double val in ratiosForCrossTalk[i])
                            sum += val;
                        writerStats.AddLine("Coverage : " + (float)(100.0 * ratiosForCrossTalk[i][i] / sum));
                        string line = "Spreading : ";
                        foreach (double val in ratiosForCrossTalk[i])
                            line += (val / sum).ToString() + ",";
                        writerStats.AddLine(line);
                    }
                    writerStats.WriteToFile();
                }
            }
        }    
    
        public static List<List<double>> ComputeMaxFlows(
                                DBOptions dbOptions,
                                Samples ProjectRatios,
                                Samples ProjectMixed,
                                List<double> RatioNormalizer,
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
            List<string> ratioNames = new List<string>();
            foreach (Sample sample in ProjectRatios)
                ratioNames.Add(sample.nameColumn);//*/

            List<double> peptideMasses = new List<double>();
            foreach (Sample sample in ProjectRatios)
                peptideMasses.Add(Peptide.ComputeMonoisotopicMass(sample.nameColumn));//*/

            //Get PSMs
            Result mixedResult = precomputedMixed;
            Result spikedResult = precomputedSpiked;

            Dictionary<int, List<List<ProductMatch>>> DicOfRatios = new Dictionary<int,List<List<ProductMatch>>>();
            Dictionary<int, List<double>>             DicOfTrackIntensity = new Dictionary<int,List<double>>();
            Dictionary<int, List<double>>             DicOfNormalizeFactor = new Dictionary<int,List<double>>();
            Dictionary<int, double>                   DicOfErrors = new Dictionary<int,double>();
            for (int nbProductsToKeep = nbProductMin; nbProductsToKeep <= nbProductMax; nbProductsToKeep++)
            {
                List<List<ProductMatch>> ratios = new List<List<ProductMatch>>();
                List<double> TrackIntensity = new List<double>();
                List<double> NormalizeFactor = new List<double>();
                BuildSinglePeptideVirtualSpectrum(spikedResult, smoothedPrecursor, nbProductsToKeep, RatioNormalizer,
                                                        ref ratios, ref TrackIntensity, ref NormalizeFactor, chargeToConsider);
                DicOfRatios.Add(nbProductsToKeep, ratios);
                DicOfTrackIntensity.Add(nbProductsToKeep, TrackIntensity);
                DicOfNormalizeFactor.Add(nbProductsToKeep, NormalizeFactor);
                DicOfErrors.Add(nbProductsToKeep, 0);
            }

            double smallestError = double.MaxValue;            
            List<List<double>> bestListOfSumOfRatio = null;
            Dictionary<Sample, List<Dictionary<int, string>>> bestDicOfResults = null;
            for (int nbProductsToKeep = nbProductMin; nbProductsToKeep <= nbProductMax; nbProductsToKeep++)
            {
                long iterError = 0;
                double cumulPercentError = 0;
                List<List<double>> listOfSumOfRatio = new List<List<double>>();
                Dictionary<Sample, List<Dictionary<int, string>>> dicOfResultsPerSample = new Dictionary<Sample, List<Dictionary<int, string>>>();
                foreach (Sample sample in ProjectMixed)
                {
                    dicOfResultsPerSample.Add(sample, new List<Dictionary<int, string>>());
                    //Dictionary<ProductSpectrum, bool> doneSpectrum = new Dictionary<ProductSpectrum, bool>();
                    //Dictionary<PeptideSpectrumMatch, ProductSpectrum> DicOfSpectrum = new Dictionary<ProductSpectrum, PeptideSpectrumMatch>();
                    Dictionary<ProductSpectrum, PeptideSpectrumMatch> DicOfSpectrum = new Dictionary<ProductSpectrum, PeptideSpectrumMatch>();

                    //Get groups of spectrum (split spectrum of different masses) into different groups
                    Dictionary<double, List<ProductSpectrum>> DicOfSpectrumMasses = new Dictionary<double, List<ProductSpectrum>>();
                    foreach (Query query in mixedResult.queries)
                    {
                        if (query.sample == sample && !DicOfSpectrum.ContainsKey(query.spectrum) && query.precursor.Charge == chargeToConsider)
                        {
                            PeptideSpectrumMatch psmToKeep = null;
                            foreach (PeptideSpectrumMatch psm in query.psms)
                                foreach (string name in ratioNames)
                                    if (name.CompareTo(psm.Peptide.Sequence) == 0)
                                        psmToKeep = psm;

                            if (psmToKeep != null)
                            {
                                bool found = false;
                                foreach (double key in DicOfSpectrumMasses.Keys)
                                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(psmToKeep.Peptide.MonoisotopicMass, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                    {
                                        found = true;
                                        DicOfSpectrumMasses[key].Add(query.spectrum);
                                    }
                                if (!found)
                                {
                                    List<ProductSpectrum> listOfSpectrum = new List<ProductSpectrum>();
                                    listOfSpectrum.Add(query.spectrum);
                                    DicOfSpectrumMasses.Add(psmToKeep.Peptide.MonoisotopicMass, listOfSpectrum);
                                }
                                DicOfSpectrum.Add(query.spectrum, psmToKeep);
                            }
                        }
                    }

                    Dictionary<int, string> dicOfResults = new Dictionary<int, string>();
                    List<double> sumOfRatio = new List<double>();
                    for (int i = 0; i < ProjectRatios.Count; i++)
                        sumOfRatio.Add(0);

                    //Cycle through groups of spectrum and evaluate metrics as if they where different samples
                    foreach (double key in DicOfSpectrumMasses.Keys)
                    {
                        List<ProductSpectrum> sortedSpectrum = DicOfSpectrumMasses[key];
                        sortedSpectrum.Sort(ProductSpectrum.AscendingRetentionTimeComparison);
                        PeptideSpectrumMatches psmMatches = new PeptideSpectrumMatches();
                        foreach (ProductSpectrum ps in sortedSpectrum)
                            psmMatches.Add(DicOfSpectrum[ps]);

                        psmMatches.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
                        Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor = psmMatches.ComputeMsMsNormalizationFactors();


                        //double LastTimeStamp = 0;
                        //double TotalElapsedTime = 0;
                        //foreach(ProductSpectrum key in sortedSpectrum)
                        foreach (PeptideSpectrumMatch psm in psmMatches)
                        {
                            //PeptideSpectrumMatch psm = DicOfSpectrum[key];

                            double overFlow = 0;
                            double underFlow = 0;
                            double percentError = 0;
                            List<List<ProductMatch>> matches = new List<List<ProductMatch>>();
                            for (int i = 0; i < peptideMasses.Count; i++)
                                if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(psm.Peptide.MonoisotopicMass, peptideMasses[i], dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                    matches.Add(DicOfRatios[nbProductsToKeep][i]);
                                else
                                    matches.Add(new List<ProductMatch>());
                            List<double> finalRatios = LaunchMaxFlowFromSpectrum(matches, ratioNames, precision, psm.Query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                                    DicOfPsmFactor[psm], mflowReturnType, ref overFlow, ref underFlow, ref percentError);

                            double normSumRatio = 0;
                            for (int i = 0; i < finalRatios.Count; i++)
                                normSumRatio += finalRatios[i] * DicOfNormalizeFactor[nbProductsToKeep][i];

                            double sumRatio = 0;
                            for (int i = 0; i < finalRatios.Count; i++)
                                sumRatio += finalRatios[i];

                            if (!double.IsNaN(normSumRatio) && normSumRatio > 0)
                            {
                                if (percentError < 0.95)
                                {
                                    List<double> avgQuantifiedRatios = new List<double>();
                                    for (int i = 0; i < ProjectRatios.Count; i++)
                                        avgQuantifiedRatios.Add((finalRatios[i] / sumRatio) * DicOfNormalizeFactor[nbProductsToKeep][i] * psm.Query.spectrum.PrecursorIntensity);

                                    string strRatios = psm.Query.spectrum.RetentionTimeInMin.ToString() + "," + psm.Query.spectrum.PrecursorIntensity;
                                    //double ElapsedTime = psm.Query.spectrum.RetentionTimeInMin - LastTimeStamp;
                                    for (int i = 0; i < avgQuantifiedRatios.Count; i++)
                                    {
                                        sumOfRatio[i] += avgQuantifiedRatios[i];
                                        strRatios += "," + avgQuantifiedRatios[i];
                                    }
                                    dicOfResults.Add(psm.Query.spectrum.ScanNumber, strRatios);
                                    //TotalElapsedTime += ElapsedTime;
                                }
                                else
                                    Console.WriteLine("Bad MaxFlow computation : " + percentError);
                            }

                            foreach (double ratio in finalRatios)
                                if (ratio == 0)
                                    percentError += 1.0 / (double)finalRatios.Count;//*/
                            cumulPercentError += percentError;
                            iterError++;
                            //LastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;
                        }
                        dicOfResultsPerSample[sample].Add(dicOfResults);
                    }
                    listOfSumOfRatio.Add(sumOfRatio);
                }
                cumulPercentError /= (double)iterError;
                if (aim4StableRatio > 0)
                {
                    foreach (List<double> listOfSum in listOfSumOfRatio)
                    {
                        double sum = 0;
                        foreach (double val in listOfSum)
                            sum += val;

                        foreach (double val in listOfSum)
                            cumulPercentError += Math.Abs(aim4StableRatio - val / sum);
                    }
                }
                
                //cumulPercentError /= (double)nbProductsToKeep;//Average error, per peak/fragment
                if (cumulPercentError < smallestError)
                {
                    nbProductsUsed = nbProductsToKeep;
                    smallestError = cumulPercentError;
                    bestListOfSumOfRatio = listOfSumOfRatio;
                    bestDicOfResults = dicOfResultsPerSample;
                }
            }
            if (bestDicOfResults != null)
            {
                Console.WriteLine("Best Number of products : " + nbProductsUsed);
                foreach (Sample sample in bestDicOfResults.Keys)
                {//Dictionary<Sample, List<Dictionary<int, string>>> 
                    vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + sample.nameColumn + "_Ratios5_Charge" + chargeToConsider + ".csv");
                    foreach (Dictionary<int, string> dic in bestDicOfResults[sample])
                    {
                        List<int> scans = new List<int>(dic.Keys);
                        scans.Sort();
                        foreach (int scan in scans)
                            writer.AddLine(dic[scan]);
                    }
                    writer.WriteToFile();
                }
            }
            return bestListOfSumOfRatio;
        }

        private static double ComputeMaxFlow(List<List<ProductMatch>> spikedMatches,
                                    List<MsMsPeak> mixedSpectrum, MassTolerance tolerance,
                                ref List<List<double>> optimalSolutions,
                                ref double percentError,
                                ref List<long> average)
        {
            //Create dictionnary of usefull peaks
            Dictionary<float, double> mixedFragDic = new Dictionary<float, double>();
            foreach (List<ProductMatch> fragmentRatio in spikedMatches)
            {
                foreach (ProductMatch match in fragmentRatio)
                {
                    if (!mixedFragDic.ContainsKey((float)match.theoMz))
                        mixedFragDic.Add((float)match.theoMz, 0);
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
                double smallestUnderError = double.MaxValue;
                double smallestOverError = overError;
                for (int i = 0; i < spikedMatches.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= iterSize;
                        virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
                        double tmpErrorMinus = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
                        double tmpErrorPlus = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);

                        if (tmpErrorPlus < overError && (tmpErrorMinus < smallestUnderError
                            || (tmpErrorMinus == smallestUnderError && tmpErrorPlus < smallestOverError)))
                        {
                            smallestOverError = tmpErrorPlus;
                            smallestUnderError = tmpErrorMinus;
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
            foreach (int localFlow in localFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            newList = new List<double>();
            foreach (long localFlow in bestLocalFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);
                                    
            //Compute average
            average.Clear();
            for (int i = 0; i < optimalSolutions[0].Count; i++)
            {
                double sum = 0.0;
                foreach (List<double> solution in optimalSolutions)
                    sum += solution[i];
                double avg = sum / (double)optimalSolutions.Count;
                average.Add((long)avg);
            }

            //Compute expected error in percentage
            virtualSpectrum = BuildVirtualSpectrum(spikedMatches, bestLocalFlows, mixedFragDic);
            double sumOfIntensities = 0;
            foreach (double val in mixedFragDic.Values)
                sumOfIntensities += val;
            percentError = underError / sumOfIntensities;

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

        private static List<double> LaunchMaxFlowFromSpectrum(List<List<ProductMatch>> ratiosToFit, List<string> ratioNames, 
                                            int precision, GraphML_List<MsMsPeak> capacity, MassTolerance tolerance, 
                                            double boostRatioForSpectrum, int returnType,//0 for max flow, 1 for best flow, 2 for average
                                            ref double overFlow, ref double underFlow, ref double errorInPercent)
        {
            List<List<double>> solutions = new List<List<double>>();
            List<long> average = new List<long>();
            List<MsMsPeak> expandedCapacity = new List<MsMsPeak>();
            double sumOfProducts = 0;
            foreach (MsMsPeak peak in capacity)
            {
                double intensity = peak.Intensity + peak.Intensity * boostRatioForSpectrum;
                expandedCapacity.Add(new MsMsPeak(peak.MZ, intensity * precision, peak.Charge));
                sumOfProducts += peak.Intensity;
            }
            List<List<ProductMatch>> tmpRatiosToFit = new List<List<ProductMatch>>();
            foreach (List<ProductMatch> list in ratiosToFit)
            {
                List<ProductMatch> pms = new List<ProductMatch>();
                foreach (ProductMatch pm in list)
                {
                    ProductMatch newPm = new ProductMatch(pm);
                    newPm.obsIntensity *= sumOfProducts;
                    pms.Add(newPm);
                }
                tmpRatiosToFit.Add(pms);
            }

            double error = ComputeMaxFlow(tmpRatiosToFit, expandedCapacity, tolerance, ref solutions, ref errorInPercent, ref average);

            double sumOfIntensities = 0;
            foreach (MsMsPeak peak in expandedCapacity)
                sumOfIntensities += peak.Intensity;

            overFlow = 0;
            underFlow = error;

            List<double> result = new List<double>();
            switch (returnType)
            {
                case 0:
                    foreach (double val in solutions[0])
                        result.Add(val / (double)precision);
                    break;
                case 1:
                    foreach (double val in solutions[1])
                        result.Add(val / (double)precision);
                    break;
                case 2:
                    foreach (double val in average)
                        result.Add(val / (double)precision);
                    break;
            }
            return result;
        }

        private static void BuildSinglePeptideVirtualSpectrum(Result precomputedResults, bool smoothPrecursor, 
                    int nbProductsToKeep, List<double> RatioNormalizer, ref List<List<ProductMatch>> FinalSpikedProducts, 
                    ref List<double> PrecursorAreas, ref List<double> Normalizor, int charge)
        {   
            DBOptions dbOptions = precomputedResults.dbOptions;
            Samples Project = precomputedResults.samples;
            
            Peptide peptide = null;
            FinalSpikedProducts.Clear();
            List<List<ProductMatch>> SpikedProducts = new List<List<ProductMatch>>();
            List<double> PeakAvgIntensities = new List<double>();
            //Get list of PSMs per sample
            double averagePrecursorArea = 0;
            double avgMsMsProductIntensity = 0;
            int nbAveragedSpectrum = 0;
            List<PeptideSpectrumMatches> listOfPSMs = new List<PeptideSpectrumMatches>();
            foreach (Sample sample in Project)
            {
                double avgPeakIntensity = 0;
                PeptideSpectrumMatches psmList = new PeptideSpectrumMatches();
                foreach (Query query in precomputedResults.queries)
                {
                    if (query.precursor.Charge == charge && query.precursor.sample == sample)
                    {
                        foreach (PeptideSpectrumMatch psm in query.psms)
                        {
                            if (sample.nameColumn.CompareTo(psm.Peptide.Sequence) == 0)// && psm.MatchingProducts > 4)//.ProbabilityScore() > 0.020)// && psm.ProbabilityScore() > bestPsmScore)
                            {
                                peptide = psm.Peptide;
                                psmList.Add(psm);
                                avgPeakIntensity += query.spectrum.PrecursorIntensity;
                                avgMsMsProductIntensity += query.spectrum.TotalIntensity;
                                nbAveragedSpectrum++;
                            }
                        }
                    }
                }

                avgPeakIntensity /= (double)psmList.Count;
                psmList.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
                listOfPSMs.Add(psmList);
                List<ProductMatch> productList = psmList.GetCombinedSpectrum(precomputedResults.dbOptions, peptide, charge);

                productList.Sort(ProductMatch.AscendingWeightComparison);

                double precursorArea = psmList.ComputePrecursorArea(smoothPrecursor);
                averagePrecursorArea += precursorArea;

                SpikedProducts.Add(productList);
                PrecursorAreas.Add(precursorArea);
                PeakAvgIntensities.Add(avgPeakIntensity);
            }
            avgMsMsProductIntensity /= nbAveragedSpectrum;
            averagePrecursorArea /= (double)Project.Count;

            Dictionary<double, int> DicOfFragmentsToKeep = new Dictionary<double,int>();
            //Get List of desired fragments, and keep masses
            for (int i = 0; i < Project.Count; i++)
            {
                for (int j = SpikedProducts[i].Count - 1; j > 0 && j > SpikedProducts[i].Count - nbProductsToKeep; j--)
                    if(!DicOfFragmentsToKeep.ContainsKey(SpikedProducts[i][j].theoMz))
                        DicOfFragmentsToKeep.Add(SpikedProducts[i][j].theoMz, 1);
                    else
                        DicOfFragmentsToKeep[SpikedProducts[i][j].theoMz] ++;
            }

            for(int i = 0; i < Project.Count; i++)
                FinalSpikedProducts.Add(listOfPSMs[i].GetCombinedSpectrum(precomputedResults.dbOptions, peptide, charge, DicOfFragmentsToKeep));
                        
            double avgPrecursors = 0;
            foreach (double avg in PeakAvgIntensities)
                avgPrecursors += avg;
            avgPrecursors /= (double)PeakAvgIntensities.Count;

            for (int i = 0; i < Project.Count; i++)
            {
                //Normalize each spectrum based on average precursor intensity
                double y = (avgPrecursors - PeakAvgIntensities[i]) / PeakAvgIntensities[i];
                for (int j = 0; j < FinalSpikedProducts[i].Count; j++)
                    FinalSpikedProducts[i][j].obsIntensity += y * FinalSpikedProducts[i][j].obsIntensity;//*/
            }
            
            List<double> listOfsumOfProducts = new List<double>();
            double avgRatioSum = 0;
            for (int i = 0; i < Project.Count; i++)
            {
                double sumOfProducts = 0;
                for (int j = 0; j < FinalSpikedProducts[i].Count; j++)
                    sumOfProducts += FinalSpikedProducts[i][j].obsIntensity;
                listOfsumOfProducts.Add(sumOfProducts);
                avgRatioSum += sumOfProducts;
            }
            avgRatioSum /= (double)Project.Count;
                        
            for (int i = 0; i < Project.Count; i++)
            {
                double spectrumWeight = 1.0;
                Normalizor.Add((spectrumWeight / listOfsumOfProducts[i]) * (averagePrecursorArea / PrecursorAreas[i]) * RatioNormalizer[i]);
                for (int j = 0; j < FinalSpikedProducts[i].Count; j++)
                    FinalSpikedProducts[i][j].obsIntensity *= (spectrumWeight / listOfsumOfProducts[i]);
            } 
        }
    }
}
