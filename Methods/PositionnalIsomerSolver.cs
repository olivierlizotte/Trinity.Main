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
                                int nbProductsMin, int nbProductsMax, bool smoothedPrecursor, int precision, int mflowReturnType, ref int nbProductsUsed)
        {
            Samples ProjectRatios = spikedResult.samples;
            Samples ProjectStable = stableResult.samples;
            List<double> tmpRatioNormalizer = new List<double>();
            foreach (Sample sample in ProjectRatios)
                tmpRatioNormalizer.Add(1.0);

            List<List<double>> tmpRatios = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectStable, tmpRatioNormalizer,
                                                            ref nbProductsUsed, mflowReturnType,
                                                            nbProductsMin, nbProductsMax, smoothedPrecursor,
                                                            precision, spikedResult, stableResult, charge);
            if (tmpRatios != null && tmpRatios.Count > 0)
            {
                List<double> RatioNormalizer = new List<double>();
                foreach (List<double> listSample in tmpRatios)
                {
                    double sum = 0;
                    int nbCoveredPrec = 0;
                    foreach (double val in listSample)
                    {
                        if (val > 0)
                            nbCoveredPrec++;

                        sum += val;
                    }
                    if (RatioNormalizer.Count == 0)
                        foreach (double val in listSample)
                            RatioNormalizer.Add(1.0);

                    double aim = 1.0 / (double)nbCoveredPrec;
                    for (int i = 0; i < listSample.Count; i++)
                        if (listSample[i] > 0)
                            RatioNormalizer[i] = aim / (listSample[i] / sum);
                }

                List<double> sampleNormalizer = new List<double>();
                foreach (List<double> listSample in tmpRatios)
                {
                    if(sampleNormalizer.Count == 0)
                        foreach (double val in listSample)
                            sampleNormalizer.Add(0);
                    for (int i = 0; i < listSample.Count; i++)
                        sampleNormalizer[i] += listSample[i] * RatioNormalizer[i];
                }
                double sumOfSamples = 0;
                for (int i = 0; i < sampleNormalizer.Count; i++)
                    sumOfSamples += sampleNormalizer[i];
                double avg = sumOfSamples / (double) sampleNormalizer.Count;
                
                for (int i = 0; i < sampleNormalizer.Count; i++)
                    RatioNormalizer[i] /= avg / sampleNormalizer[i];
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

                List<double> RatioNormalizer = null;
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
                    RatioNormalizer = new List<double>();
                    for (int i = 0; i < ProjectRatios.Count; i++)
                        RatioNormalizer.Add(1.0);
                }

                //Solve positional isomers in the mixed spectrums
                List<List<double>> ratios = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectMixed, RatioNormalizer, ref nbProductsUsed, mflowReturnType,
                                                            nbProductsMin, nbProductsMax, smoothedPrecursor, precision, spikedResult, mixedResult, charge);

                if (ratios == null || ratios.Count == 0)
                    dbOptions.ConSole.WriteLine("No precursor found with charge " + charge);
                else
                {
                    //Export results to a csv file
                    vsCSVWriter writerCumul = new vsCSVWriter(dbOptions.OutputFolder + "CumulRatios9_charge" + charge + ".csv");
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
                    //writerStats.AddLine("aim (for stable ratios) : " + aim);
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
            //Export results of searches
            spikedResult.Export(0.05);
            mixedResult.Export(0.05);
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

            Dictionary<int, Dictionary<double, List<List<ProductMatch>>>>   DicOfRatios = new Dictionary<int, Dictionary<double, List<List<ProductMatch>>>>();
            Dictionary<int, Dictionary<double, List<double>>>               DicOfTrackIntensity = new Dictionary<int, Dictionary<double, List<double>>>();
            Dictionary<int, Dictionary<double, List<double>>>               DicOfNormalizeFactor = new Dictionary<int, Dictionary<double, List<double>>>();
            Dictionary<int, double>                                         DicOfErrors = new Dictionary<int,double>();
            Dictionary<int, List<double>>                                   DicOfFragmentMz = new Dictionary<int, List<double>>();

            for (int nbProductsToKeep = nbProductMin; nbProductsToKeep <= nbProductMax; nbProductsToKeep++)
            {
                Dictionary<double, List<List<ProductMatch>>> ratios = new Dictionary<double, List<List<ProductMatch>>>();
                Dictionary<double, List<double>> TrackIntensity = new Dictionary<double, List<double>>();
                Dictionary<double, List<double>> NormalizeFactor = new Dictionary<double, List<double>>();
                List<double> ListFragmentMz = new List<double>();

                BuildSinglePeptideVirtualSpectrum(spikedResult, smoothedPrecursor, nbProductsToKeep, RatioNormalizer,
                                                        ref ratios, ref TrackIntensity, ref NormalizeFactor, ref ListFragmentMz, chargeToConsider);

                DicOfRatios.Add(nbProductsToKeep, ratios);
                DicOfTrackIntensity.Add(nbProductsToKeep, TrackIntensity);
                DicOfNormalizeFactor.Add(nbProductsToKeep, NormalizeFactor);
                DicOfFragmentMz.Add(nbProductsToKeep, ListFragmentMz);
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
                    List<double> listOfSumOfRatioPerPrecursor = new List<double>();
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
                        List<double> sumOfRatio = new List<double>();
                        for (int i = 0; i < ProjectRatios.Count; i++)
                            sumOfRatio.Add(0);

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


                        Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor = psmMatches.ComputeMsMsNormalizationFactors();
                        
                        Dictionary<int, string> dicOfResults = new Dictionary<int, string>();
                        
                        //double precursorArea = psmMatches.ComputePrecursorArea(smoothedPrecursor);
                        double precursorArea = 0.0;
                        double LastIntensityPerUnitOrTime = 0;
                        double LastTimeStamp = 0;
                        //double TotalElapsedTime = 0;
                        //foreach(ProductSpectrum key in sortedSpectrum)
                        foreach (PeptideSpectrumMatch psm in psmMatches)
                        {
                            double overFlow = 0;
                            double underFlow = 0;
                            double percentError = 0;
                            List<List<ProductMatch>> matches = new List<List<ProductMatch>>();
                            for (int i = 0; i < peptideMasses.Count; i++)
                                if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(psm.Peptide.MonoisotopicMass, peptideMasses[i], dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                    matches.Add(DicOfRatios[nbProductsToKeep][key][i]);
                                else
                                    matches.Add(new List<ProductMatch>());
                            List<double> finalRatios = LaunchMaxFlowFromSpectrum(matches, ratioNames, precision, psm.Query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                                    DicOfPsmFactor[psm], mflowReturnType, ref overFlow, ref underFlow, ref percentError, dbOptions.ConSole);

                            double normSumRatio = 0;
                            for (int i = 0; i < finalRatios.Count; i++)
                                normSumRatio += finalRatios[i] * DicOfNormalizeFactor[nbProductsToKeep][key][i];

                            double sumRatio = 1.0;
                            //for (int i = 0; i < finalRatios.Count; i++)
                            //    sumRatio += finalRatios[i];

                            double ElapsedTime = (psm.Query.spectrum.RetentionTimeInMin - LastTimeStamp) * 60.0 * 1000.0;
                            double IntensityPerUnitOfTime = psm.Query.spectrum.PrecursorIntensity / psm.Query.spectrum.Ms1InjectionTime;// * ElapsedTime;
                                   
                            double localArea = (IntensityPerUnitOfTime + LastIntensityPerUnitOrTime) * 0.5 * ElapsedTime;                            

                            if (LastTimeStamp > 0)
                            {
                                //if (IntensityPerUnitOfTime > 1.75 * LastIntensityPerUnitOrTime)
                                //    Console.WriteLine("oops?");
                                precursorArea += localArea;
                                if (percentError < 0.25 && !double.IsNaN(normSumRatio) && normSumRatio > 0)//0.95
                                {
                                    List<double> avgQuantifiedRatios = new List<double>();

                                    for (int i = 0; i < ProjectRatios.Count; i++)
                                    {
                                        double avgRatio = (finalRatios[i] / sumRatio) * DicOfNormalizeFactor[nbProductsToKeep][key][i] * localArea;
                                        if (double.IsNaN(avgRatio))
                                            dbOptions.ConSole.WriteLine("Oops, NaN in ratios");
                                        avgQuantifiedRatios.Add(avgRatio);
                                    }

                                    string strRatios = psm.Query.spectrum.RetentionTimeInMin.ToString() + "," + localArea;

                                    for (int i = 0; i < avgQuantifiedRatios.Count; i++)
                                    {
                                        sumOfRatio[i] += avgQuantifiedRatios[i];
                                        strRatios += "," + avgQuantifiedRatios[i];
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

                                foreach (double ratio in finalRatios)
                                    if (ratio == 0)
                                        percentError += 1.0 / (double)finalRatios.Count;//*/
                                cumulPercentError += percentError;
                                iterError++;
                            }

                            LastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;     
                            LastIntensityPerUnitOrTime = IntensityPerUnitOfTime;
                        }

                        dicOfResultsPerSample[sample].Add(dicOfResults);

                        //Use SumOfRatio as the fraction for the entire Precursor area
                        /*double sumArea = 0;
                        foreach (double val in sumOfRatio)
                            sumArea += val;
                        
                        List<double> relativeSumOfRatio = new List<double>();
                        for (int i = 0; i < sumOfRatio.Count; i++)
                            sumOfRatio[i] = (sumOfRatio[i] / sumArea) * precursorArea;
                        //*/
                        if (listOfSumOfRatioPerPrecursor.Count == 0)
                            foreach (double area in sumOfRatio)
                                listOfSumOfRatioPerPrecursor.Add(area);
                        else
                            for (int i = 0; i < sumOfRatio.Count; i++)
                                listOfSumOfRatioPerPrecursor[i] += sumOfRatio[i];
                    }//precursor

                    listOfSumOfRatio.Add(listOfSumOfRatioPerPrecursor);//Per Precursor

                }//sample

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
                dbOptions.ConSole.WriteLine("Best Number of products : " + nbProductsUsed);
                foreach (Sample sample in bestDicOfResults.Keys)
                {//Dictionary<Sample, List<Dictionary<int, string>>> 
                    vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + sample.nameColumn + "_Ratios5_Charge" + chargeToConsider + ".csv");
                    string title = "Retention Time,Precursor Intensity";
                    for(int i = 0; i < ProjectRatios.Count; i++)
                        title += "," + ProjectRatios[i].nameColumn;
                    
                    foreach (double mz in DicOfFragmentMz[nbProductsUsed])
                        title += "," + mz;
                    writer.AddLine(title);
                    foreach (Dictionary<int, string> dic in bestDicOfResults[sample])
                    {
                        List<int> scans = new List<int>(dic.Keys);
                        scans.Sort();
                        foreach (int scan in scans)
                            writer.AddLine(dic[scan]);
                    }
                    writer.WriteToFile();
                }
                foreach (List<double> sampleRatio in bestListOfSumOfRatio)
                    for (int i = 0; i < sampleRatio.Count; i++)
                        sampleRatio[i] *= RatioNormalizer[i];
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
                                            ref double overFlow, ref double underFlow, ref double errorInPercent, IConSol ConSole)
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
            return result;
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
                    int nbProductsToKeep, List<double> RatioNormalizer, ref Dictionary<double, List<List<ProductMatch>>> FinalSpikedProducts, 
                    ref Dictionary<double, List<double>> PrecursorAreas, ref Dictionary<double, List<double>> Normalizor, ref List<double> FragmentMz, int charge)
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
                    FinalSpikedProducts.Add(foundKey, new List<List<ProductMatch>>());
                    PrecursorAreas.Add(foundKey, new List<double>());
                    Normalizor.Add(foundKey, new List<double>());
                }
            }
            
            Dictionary<double, int> AllFragments = new Dictionary<double,int>();
            vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + "FragmentsUsed_" + nbProductsToKeep + "Products.csv");
            writer.AddLine("Precursor Mass, Peptide Sequence,Fragment,Pos,Charge,Mz,Intensity");
            foreach (double foundKey in FinalSpikedProducts.Keys)
            {
                Peptide peptide = null;
                List<List<ProductMatch>> SpikedProducts = new List<List<ProductMatch>>();
                List<double> PeakAvgIntensities = new List<double>();
                //Get list of PSMs per sample
                double averagePrecursorArea = 0;
                double avgMsMsProductIntensity = 0;
                int nbAveragedSpectrum = 0;
                int nbAverageArea = 0;
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
                                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(psm.Peptide.MonoisotopicMass, foundKey, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
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
                    }

                    if (psmList.Count > 0)
                    {
                        avgPeakIntensity /= (double)psmList.Count;
                        psmList.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
                        listOfPSMs.Add(psmList);
                        List<ProductMatch> productList = psmList.GetCombinedSpectrum(precomputedResults.dbOptions, peptide, charge);

                        productList.Sort(ProductMatch.AscendingWeightComparison);

                        double precursorArea = psmList.ComputePrecursorArea(smoothPrecursor);
                        averagePrecursorArea += precursorArea;
                        nbAverageArea++;

                        SpikedProducts.Add(productList);
                        PrecursorAreas[foundKey].Add(precursorArea);
                        PeakAvgIntensities.Add(avgPeakIntensity);
                    }
                    else
                    {
                        listOfPSMs.Add(psmList);
                        SpikedProducts.Add(new List<ProductMatch>());
                        PrecursorAreas[foundKey].Add(0);
                        PeakAvgIntensities.Add(0);
                    }
                }
                avgMsMsProductIntensity /= nbAveragedSpectrum;
                averagePrecursorArea /= (double)nbAverageArea;

                Dictionary<double, int> DicOfFragmentsToKeep = new Dictionary<double, int>();
                //Get List of desired fragments, and keep masses
                for (int i = 0; i < Project.Count; i++)
                {
                    for (int j = SpikedProducts[i].Count - 1; j > 0 && j >= SpikedProducts[i].Count - nbProductsToKeep; j--)
                        if (!DicOfFragmentsToKeep.ContainsKey(SpikedProducts[i][j].theoMz))
                            DicOfFragmentsToKeep.Add(SpikedProducts[i][j].theoMz, 1);
                        else
                            DicOfFragmentsToKeep[SpikedProducts[i][j].theoMz]++;
                }

                for (int i = 0; i < Project.Count; i++)
                {
                    List<ProductMatch> list = listOfPSMs[i].GetCombinedSpectrum(precomputedResults.dbOptions, peptide, charge, DicOfFragmentsToKeep);
                                        
                    foreach (ProductMatch match in list)
                        writer.AddLine(foundKey + "," + Project[i].nameColumn + "," + match.fragment + "," + match.fragmentPos + "," + match.charge + "," + match.theoMz + "," + match.obsIntensity);

                    FinalSpikedProducts[foundKey].Add(list);
                }

                foreach(double key in DicOfFragmentsToKeep.Keys)
                    if(!AllFragments.ContainsKey(key))
                        AllFragments.Add(key, 1);
                    else
                        AllFragments[key]++;

                double avgPrecursors = 0;
                int nbInt = 0;
                foreach (double avg in PeakAvgIntensities)
                {
                    if (avg > 0)
                        nbInt++;
                    avgPrecursors += avg;
                }
                avgPrecursors /= (double)nbInt;

                for (int i = 0; i < Project.Count; i++)
                {
                    //Normalize each spectrum based on average precursor intensity
                    if (PeakAvgIntensities[i] > 0)
                    {
                        double y = (avgPrecursors - PeakAvgIntensities[i]) / PeakAvgIntensities[i];
                        for (int j = 0; j < FinalSpikedProducts[foundKey][i].Count; j++)
                            FinalSpikedProducts[foundKey][i][j].obsIntensity += y * FinalSpikedProducts[foundKey][i][j].obsIntensity;//*/
                    }
                }

                List<double> listOfsumOfProducts = new List<double>();
                double avgRatioSum = 0;
                nbInt = 0;
                for (int i = 0; i < Project.Count; i++)
                {
                    double sumOfProducts = 0;
                    for (int j = 0; j < FinalSpikedProducts[foundKey][i].Count; j++)
                        sumOfProducts += FinalSpikedProducts[foundKey][i][j].obsIntensity;
                    listOfsumOfProducts.Add(sumOfProducts);

                    avgRatioSum += sumOfProducts; 

                    if (sumOfProducts > 0)
                        nbInt++;
                }
                avgRatioSum /= (double)nbInt;

                for (int i = 0; i < Project.Count; i++)
                {
                    if (listOfsumOfProducts[i] > 0)
                    {
                        double lossOfPrecIntensity = Math.Log(averagePrecursorArea, 2) / Math.Log(PrecursorAreas[foundKey][i], 2);
                        double lossofFragmentIntensity = Math.Pow(avgRatioSum, 2) / Math.Pow(listOfsumOfProducts[i], 2);
                        //Normalizor[foundKey].Add((spectrumWeight / listOfsumOfProducts[i]) * (averagePrecursorArea / PrecursorAreas[foundKey][i]) * RatioNormalizer[i]);
                        //Normalizor[foundKey].Add(lossOfPrecIntensity * lossofFragmentIntensity);//(averagePrecursorArea / PrecursorAreas[foundKey][i]) * RatioNormalizer[i]);
                        Normalizor[foundKey].Add(lossofFragmentIntensity * lossOfPrecIntensity);//(averagePrecursorArea / PrecursorAreas[foundKey][i]) * RatioNormalizer[i]);
                        for (int j = 0; j < FinalSpikedProducts[foundKey][i].Count; j++)
                            FinalSpikedProducts[foundKey][i][j].obsIntensity /= listOfsumOfProducts[i];
                    }
                    else
                        Normalizor[foundKey].Add(1.0);
                }
            }

            foreach(double key in AllFragments.Keys)
                FragmentMz.Add(key);

            writer.WriteToFile();
        }
    }
}
