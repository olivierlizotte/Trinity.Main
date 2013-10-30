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
    public class IsomericPeptideQuantifier2
    {
        public static DBOptions GetDBOptions(bool loadFromRaw, bool onlyYions)
        {
            string outputDir = @"C:\_IRIC\DATA\Test\testNB\Iso2\";
            string fastaFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\FEB13_2013\MSMS files\Peptide.fasta";
            DBOptions dbOptions = new DBOptions(fastaFile);
            dbOptions.precursorMassTolerance = new MassTolerance(8/*8*//*8withoutisotopes*/, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(20/*8*//*8withoutisotopes*/, MassToleranceUnits.ppm);
            //dbOptions.productMassTolerance = new MassTolerance(0.05/*0.034*//*without isotopes*/, MassToleranceUnits.Da);//0.034 is a 60 000 resolution over 2000 range in mz
            dbOptions.MaximumPeptideMass = 200000;
            dbOptions.OutputFolder = outputDir;
            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            dbOptions.DigestionEnzyme = proteases["no enzyme"];//trypsin (no proline rule)"];
            dbOptions.NoEnzymeSearch = false;// true;
            dbOptions.DecoyFusion = false;
            dbOptions.MaximumNumberOfFragmentsPerSpectrum = 400;
            //dbOptions.protease = proteases["trypsin (no proline rule)"];
            dbOptions.ToleratedMissedCleavages = 200;// 2;
            dbOptions.MinimumPeptideLength = 5;
            dbOptions.MaximumPeptideLength = 300;

            GraphML_List<Modification> fixMods = new GraphML_List<Modification>();
            //fixMods.Add(ModificationDictionary.Instance["propionylation of K"]);
            dbOptions.fixedModifications = fixMods;

            GraphML_List<Modification> varMods = new GraphML_List<Modification>();
            varMods.Add(ModificationDictionary.Instance["acetylation of K"]);
            varMods.Add(ModificationDictionary.Instance["propionylation of K"]);
            dbOptions.maximumVariableModificationIsoforms = 1024;
            dbOptions.variableModifications = varMods;

            dbOptions.addFragmentLoss = false;// true;
            dbOptions.addFragmentMods = false;// true;
            dbOptions.fragments = new Fragments();
            if (!onlyYions)
            {
                dbOptions.fragments.Add(new FragmentA());
                dbOptions.fragments.Add(new FragmentB());
                dbOptions.fragments.Add(new FragmentC());
                dbOptions.fragments.Add(new FragmentX());
                dbOptions.fragments.Add(new FragmentZ());
            }
            dbOptions.fragments.Add(new FragmentY());

            dbOptions.SaveMS1Peaks = true;
            dbOptions.SaveMSMSPeaks = true;
            dbOptions.LoadSpectraIfFound = !loadFromRaw;

            dbOptions.NbPSMToKeep = 100;
            return dbOptions;
        }

        public static void MaxFlowThis()//string projectSingleInjections, string projectMixed)
        {
            //0.0732351971695529 : Return type (1) Precision (1000000) nbProductsToKeep (5)
            //0.0289455865219887 : Return type (1) Precision (100000) nbProductsToKeep (13)
            Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            Samples ProjectMixed =  new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Varied_19Oct.csv", 0);
            Samples ProjectStable = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_MonoAce_19Oct.csv", 0);

            //0.305567296686846 : Return type (0) Precision (100000) nbProductsToKeep (5)
            //Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_DiAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            //Samples ProjectMixed =  new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_DiAce_Varied_19Oct.csv", 0);
            //Samples ProjectStable = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_DiAce_19Oct.csv", 0);

            //0.222529420854448 : Return type (2) Precision (100000) nbProductsToKeep (19)
            //Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_TriAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            //Samples ProjectMixed =  new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_TriAce_Varied_19Oct.csv", 0);
            //Samples ProjectStable = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_TriAce_19Oct.csv", 0);

            DBOptions dbOptions = GetDBOptions(false, false);
            int mflowReturnType = 0;
            int nbProductsMin = 4;// 4;
            int nbProductsMax = 14;// 14;
            bool smoothedPrecursor = false;
            int precision = 1000;
            int nbProductsUsed = 0;

            double aim = 1.0 / (double)ProjectRatios.Count;

            Result spikedResult = null;
            {
                Propheus propheus = new Propheus(dbOptions, ProjectRatios);
                propheus.Preload(false, false);
                propheus.PrepareQueries();

                spikedResult = propheus.SearchLatestVersion(propheus.AllQueries, false);
            }
            List<double> RatioNormalizer = new List<double>();
            foreach (Sample sample in ProjectRatios)
                RatioNormalizer.Add(1.0);
            Result stableResult = null;
            {
                Propheus propheus = new Propheus(dbOptions, ProjectStable);
                propheus.Preload(false, false);
                propheus.PrepareQueries();

                stableResult = propheus.SearchLatestVersion(propheus.AllQueries, false);

                double smallestNormalization = double.MaxValue;
                //for (int i = nbProductsMin; i <= nbProductsMax; i++)
                {
                    int tmpNbProductsUsed = 0;
                    List<double> tmpRatioNormalizer = new List<double>();
                    foreach (Sample sample in ProjectRatios)
                        tmpRatioNormalizer.Add(1.0);
                    List<List<double>> tmpRatios = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectStable, tmpRatioNormalizer,
                                                                    ref tmpNbProductsUsed, mflowReturnType, 
                                                                    nbProductsMin, nbProductsMax, smoothedPrecursor, 
                                                                    precision, spikedResult, stableResult, aim);

                    tmpRatioNormalizer.Clear();
                    double sum = 0;
                    foreach (double val in tmpRatios[0])
                        sum += val;

                    double sumOfCorrection = 0;
                    foreach (double val in tmpRatios[0])
                        sumOfCorrection += Math.Abs(aim - val / sum);
                    if (sumOfCorrection < smallestNormalization)
                    {
                        nbProductsUsed = tmpNbProductsUsed;
                        smallestNormalization = sumOfCorrection;
                        RatioNormalizer.Clear();
                        foreach (double val in tmpRatios[0])
                            RatioNormalizer.Add(aim / (val / sum));
                    }
                }
            }

            nbProductsMin = nbProductsUsed;
            nbProductsMax = nbProductsUsed;
            List<List<double>> ratios = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectMixed, RatioNormalizer, ref nbProductsUsed,
                                                        mflowReturnType, 
                                                        nbProductsMin, nbProductsMax, smoothedPrecursor, precision, spikedResult);

            vsCSVWriter writerCumul = new vsCSVWriter(dbOptions.OutputFolder + "CumulRatios7.csv");
            for (int i = 0; i < ProjectMixed.Count; i++)
            {
                string lineCumulRatio = ProjectMixed[i].sSDF;
                for (int k = 0; k < ratios[i].Count; k++)
                    lineCumulRatio += "," + ratios[i][k];// *RatioNormalizer[k];
                writerCumul.AddLine(lineCumulRatio);
            }
            //writerCumul.AddLine(strNorm);
            writerCumul.WriteToFile();

            vsCSVWriter writerStats = new vsCSVWriter(dbOptions.OutputFolder + "Stats.txt");
            writerStats.AddLine("ProjectRatios : " + ProjectRatios.FileName);
            writerStats.AddLine("ProjectStable : " + ProjectStable.FileName);
            writerStats.AddLine("ProjectMixed : " + ProjectMixed.FileName);
            writerStats.AddLine("Nb Of products considered : " + nbProductsUsed);
            writerStats.AddLine("Precision : " + precision);
            writerStats.AddLine("mflowReturnType : " + mflowReturnType);
            writerStats.AddLine("smoothedPrecursor : " + smoothedPrecursor);
            writerStats.AddLine("aim (for stable ratios) : " + aim);
            for (int i = 0; i < RatioNormalizer.Count; i++)
                writerStats.AddLine("Normalization for " + vsCSV.GetFileName_NoExtension(ProjectRatios[i].sSDF) + " : " + RatioNormalizer[i]);

            writerStats.AddLine(" -- Cross Talk between Spiked Samples -- ");
            List<List<double>> ratiosForCrossTalk = ComputeMaxFlows(dbOptions, ProjectRatios, ProjectRatios, RatioNormalizer, ref nbProductsUsed, mflowReturnType,
                                                                    nbProductsMin, nbProductsMax, smoothedPrecursor, precision, spikedResult);

            for (int i = 0; i < ratiosForCrossTalk.Count; i++)
            {
                writerStats.AddLine(vsCSV.GetFileName_NoExtension(ProjectRatios[i].sSDF));
                double sum = 0;
                foreach (double val in ratiosForCrossTalk[i])
                    sum += val;
                writerStats.AddLine("Coverage : " + (float)(100.0 * ratiosForCrossTalk[i][i] / sum));
                string line = "Spreading : ";
                foreach (double val in ratiosForCrossTalk[i])
                    line += (val/sum).ToString() + ",";
                writerStats.AddLine(line);
            }
            writerStats.WriteToFile();
        }    
    
        public static List<List<double>> ComputeMaxFlows(
                                DBOptions dbOptions,
                                Samples ProjectRatios,
                                Samples ProjectMixed,
                                List<double> RatioNormalizer,
                                ref int nbProductsUsed,
                                int mflowReturnType = 0,
                                int nbProductMin = 5,// 5,
                                int nbProductMax = 15,// 5,
                                bool smoothedPrecursor = false,
                                int precision = 100000,
                                Result precomputedSpiked = null,
                                Result precomputedMixed  = null,
                                double aim4StableRatio = -1)
        {
            string baseSeq = "";
                        
            List<string> ratioNames = new List<string>();
            foreach (Sample sample in ProjectRatios)
                ratioNames.Add(sample.nameColumn);//*/

            //Get PSMs
            Result mixedResult = precomputedMixed;
            if (precomputedMixed == null)
            {
                Propheus propheus = new Propheus(dbOptions, ProjectMixed);

                propheus.Preload(false, false);
                propheus.PrepareQueries();

                mixedResult = propheus.SearchLatestVersion(propheus.AllQueries, false);
            }

            Result spikedResult = precomputedSpiked;
            if (precomputedSpiked == null)
            {
                Propheus propheusSpike = new Propheus(dbOptions, ProjectRatios);

                propheusSpike.Preload(false, false);
                propheusSpike.PrepareQueries();

                spikedResult = propheusSpike.SearchLatestVersion(propheusSpike.AllQueries, false);
            }

            Dictionary<int, List<List<ProductMatch>>> DicOfRatios = new Dictionary<int,List<List<ProductMatch>>>();
            Dictionary<int, List<double>>             DicOfTrackIntensity = new Dictionary<int,List<double>>();
            Dictionary<int, List<double>>             DicOfNormalizeFactor = new Dictionary<int,List<double>>();
            Dictionary<int, List<List<double>>>       DicOfComputedRatios = new Dictionary<int,List<List<double>>>();
            Dictionary<int, double>                   DicOfErrors = new Dictionary<int,double>();
            for (int nbProductsToKeep = nbProductMin; nbProductsToKeep <= nbProductMax; nbProductsToKeep++)
            {
                List<List<ProductMatch>> ratios = new List<List<ProductMatch>>();
                List<double> TrackIntensity = new List<double>();
                List<double> NormalizeFactor = new List<double>();
                CreateVirtualSpectrumFromSpikedPeptides(dbOptions, ProjectRatios, ref baseSeq, false, nbProductsToKeep, RatioNormalizer,
                                                        ref ratios, ref TrackIntensity, ref NormalizeFactor, smoothedPrecursor, spikedResult);
                DicOfRatios.Add(nbProductsToKeep, ratios);
                DicOfTrackIntensity.Add(nbProductsToKeep, TrackIntensity);
                DicOfNormalizeFactor.Add(nbProductsToKeep, NormalizeFactor);
                DicOfComputedRatios.Add(nbProductsToKeep, new List<List<double>>());
                DicOfErrors.Add(nbProductsToKeep, 0);
            }

            double smallestError = double.MaxValue;            
            List<List<double>> bestListOfSumOfRatio = null;
            Dictionary<Sample, Dictionary<int, string>> bestDicOfResults = null;
            for (int nbProductsToKeep = nbProductMin; nbProductsToKeep <= nbProductMax; nbProductsToKeep++)
            {
                long iterError = 0;
                double cumulPercentError = 0;
                List<List<double>> listOfSumOfRatio = new List<List<double>>();
                Dictionary<Sample, Dictionary<int, string>> dicOfResultsPerSample = new Dictionary<Sample, Dictionary<int, string>>();
                foreach(Sample sample in ProjectMixed)
                {
                    //Console.WriteLine("Sample " + sample.sSDF);

                    Dictionary<ProductSpectrum, bool> doneSpectrum = new Dictionary<ProductSpectrum, bool>();
                    Dictionary<int, string> dicOfResults = new Dictionary<int, string>();
                    List<double> sumOfRatio = new List<double>();
                    for (int i = 0; i < ProjectRatios.Count; i++ )
                        sumOfRatio.Add(0);

                    Dictionary<ProductSpectrum, PeptideSpectrumMatch> DicOfSpectrum = new Dictionary<ProductSpectrum, PeptideSpectrumMatch>();
                    foreach (Query query in mixedResult.queries)
                    {
                        if (query.sample == sample)
                        {
                            foreach (PeptideSpectrumMatch psm in query.psms)
                            {
                                bool keep = false;
                                foreach (string name in ratioNames)
                                    if (name.CompareTo(psm.Peptide.Sequence) == 0)
                                        keep = true;

                                if (keep && !DicOfSpectrum.ContainsKey(query.spectrum) && query.precursor.Charge == 2)// && psm.MatchingProducts > 4)
                                    DicOfSpectrum.Add(query.spectrum, psm);
                            }
                        }
                    }
                    List<ProductSpectrum> sortedSpectrum = new List<ProductSpectrum>(DicOfSpectrum.Keys);
                    sortedSpectrum.Sort(ProductSpectrum.AscendingRetentionTimeComparison);
                    //double LastTimeStamp = 0;
                    //double TotalElapsedTime = 0;
                    foreach(ProductSpectrum key in sortedSpectrum)
                    {
                        PeptideSpectrumMatch psm = DicOfSpectrum[key];

                        double overFlow = 0;
                        double underFlow = 0;
                        double percentError = 0;
                        List<double> finalRatios = MaxFlowFromSpectrum(DicOfRatios[nbProductsToKeep], ratioNames, precision, psm.Query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                                mflowReturnType, ref overFlow, ref underFlow, ref percentError);
                                                
                        double normSumRatio = 0;
                        for (int i = 0; i < finalRatios.Count; i++)
                            normSumRatio += finalRatios[i] * DicOfNormalizeFactor[nbProductsToKeep][i];

                        double sumRatio = 0;
                        for (int i = 0; i < finalRatios.Count; i++)
                            sumRatio += finalRatios[i];
                        
                        if (normSumRatio != double.NaN && normSumRatio > 0)// && nbProductsToKeep == 5)
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
                    dicOfResultsPerSample.Add(sample, dicOfResults);
                    listOfSumOfRatio.Add(sumOfRatio);
                }
                cumulPercentError /= (double)iterError;
                if (aim4StableRatio > 0)
                {
                    List<double> tmpRatioNormalizer = new List<double>();
                    double sum = 0;
                    foreach (double val in listOfSumOfRatio[0])
                        sum += val;
                    double aim = 1.0 / (double)listOfSumOfRatio[0].Count;

                    foreach (double val in listOfSumOfRatio[0])
                        cumulPercentError += Math.Abs(aim - val / sum);
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
                {
                    vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + sample.nameColumn + "_Ratios5.csv");
                    List<int> scans = new List<int>(bestDicOfResults[sample].Keys);
                    scans.Sort();
                    foreach (int scan in scans)
                        writer.AddLine(bestDicOfResults[sample][scan]);
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
            double overError = ComputeOverflow(virtualSpectrum, mixedFragDic);
            double underError = ComputeUnderflow(virtualSpectrum, mixedFragDic);
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
                        double tmpErrorMinus = ComputeUnderflow(virtualSpectrum, mixedFragDic);
                        double tmpErrorPlus = ComputeOverflow(virtualSpectrum, mixedFragDic);

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
                overError = ComputeOverflow(virtualSpectrum, mixedFragDic);
                underError = ComputeUnderflow(virtualSpectrum, mixedFragDic);
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
            
            return ComputeUnderflow(virtualSpectrum, mixedFragDic);
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

        private static double ComputeOverflow(Dictionary<float, double> virtualFragments, Dictionary<float, double> realFragments)
        {
            double overError = 0;
            foreach (float key in virtualFragments.Keys)
            {
                double sum = virtualFragments[key] - realFragments[key];
                if (sum > 0)
                    overError += sum;
            }
            return overError;
        }

        private static double ComputeUnderflow(Dictionary<float, double> virtualFragments, Dictionary<float, double> realFragments)
        {

            double underError = 0;
            foreach (float key in virtualFragments.Keys)
            {
                double sum = realFragments[key] - virtualFragments[key];
                if (sum > 0)
                    underError += sum;
            }
            return underError;
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

        public static List<double> MaxFlowFromSpectrum(List<List<ProductMatch>> ratiosToFit, List<string> ratioNames, 
                                            int precision, GraphML_List<MsMsPeak> capacity, MassTolerance tolerance, 
                                            int returnType,//0 for max flow, 1 for best flow, 2 for average
                                            ref double overFlow, ref double underFlow, ref double errorInPercent)
        {
            List<List<double>> solutions = new List<List<double>>();
            List<long> average = new List<long>();
            List<MsMsPeak> expandedCapacity = new List<MsMsPeak>();
            double sumOfProducts = 0;
            foreach (MsMsPeak peak in capacity)
            {
                expandedCapacity.Add(new MsMsPeak(peak.MZ, peak.Intensity * precision, peak.Charge));
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

        public static void CreateVirtualSpectrumFromSpikedPeptides(DBOptions dbOptions, Samples Project, ref string BaseSequence, bool loadFromRaw, int nbProductsToKeep,
                    List<double> RatioNormalizer,
                    ref List<List<ProductMatch>> FinalSpikedProducts, ref List<double> PrecursorAreas, ref List<double> Normalizor, bool smoothPrecursor, Result precomputed)
        {   
            int charge = 2;
            Peptide peptide = null;
            FinalSpikedProducts.Clear();
            List<List<ProductMatch>> SpikedProducts = new List<List<ProductMatch>>();
            List<double> PeakAvgIntensities = new List<double>();
            //Get list of PSMs per sample
            double averagePrecursorArea = 0;
            double avgMsMsProductIntensity = 0;
            int nbAveragedSpectrum = 0;
            List<List<PeptideSpectrumMatch>> listOfPSMs = new List<List<PeptideSpectrumMatch>>();
            foreach (Sample sample in Project)
            {
                double avgPeakIntensity = 0;
                PeptideSpectrumMatches psmList = new PeptideSpectrumMatches();
                foreach (Query query in precomputed.queries)
                {
                    if (query.precursor.Charge == charge && query.precursor.sample == sample)
                    {
                        foreach (PeptideSpectrumMatch psm in query.psms)
                        {
                            if (sample.nameColumn.CompareTo(psm.Peptide.Sequence) == 0)// && psm.MatchingProducts > 4)//.ProbabilityScore() > 0.020)// && psm.ProbabilityScore() > bestPsmScore)
                            {
                                peptide = psm.Peptide;
                                BaseSequence = peptide.BaseSequence;
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
                List<ProductMatch> productList = GetCombinedSpectrum(precomputed.dbOptions, psmList, peptide, charge);

                productList.Sort(ProductMatch.AscendingWeightComparison);
                //if (productList.Count > nbProductsToKeep)
                //    productList.RemoveRange(0, productList.Count - nbProductsToKeep);//*/

                double precursorArea = psmList.ComputePrecursorArea(smoothPrecursor);
                averagePrecursorArea += precursorArea;

                SpikedProducts.Add(productList);
                PrecursorAreas.Add(precursorArea);// / (double)nbPrec);
                PeakAvgIntensities.Add(avgPeakIntensity);
                //tmp.ExportFragmentIntensitiesForAllPSM(psmList, peptide, charge, tmp.dbOptions.OutputFolder + sample.nameColumn + ".csv");
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
                FinalSpikedProducts.Add(GetCombinedSpectrum(precomputed.dbOptions, listOfPSMs[i], peptide, charge, DicOfFragmentsToKeep));
                        
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
            /*
            double avgSumProducts = 0;
            foreach (double sum in listOfsumOfProducts)
                avgSumProducts += sum;
            avgSumProducts /= (double)listOfsumOfProducts.Count;//*/
                        
            for (int i = 0; i < Project.Count; i++)
            {
                double spectrumWeight = 1.0;// avgMsMsProductIntensity;// 10000000.0;
                Normalizor.Add((spectrumWeight / listOfsumOfProducts[i]) * (averagePrecursorArea / PrecursorAreas[i]) * RatioNormalizer[i]);
                for (int j = 0; j < FinalSpikedProducts[i].Count; j++)
                    FinalSpikedProducts[i][j].obsIntensity *= (spectrumWeight / listOfsumOfProducts[i]);
            } 
        }

        public static Dictionary<PeptideSpectrumMatch, double> ComputeMsMsNormalizationFactors(List<PeptideSpectrumMatch> psms)
        {
            Dictionary<PeptideSpectrumMatch, double> fragRatio = new Dictionary<PeptideSpectrumMatch, double>();
            double lastIntensity = 0;
            foreach (PeptideSpectrumMatch psm in psms)
            {
                if (psm.Query.spectrum.PrecursorIntensity > 0)
                {
                    double intensityFactor = 0;
                    if (psm.Query.spectrum.InjectionTime >= 119.999997317791)
                        lastIntensity = psm.Query.spectrum.PrecursorIntensity;
                    else
                    {
                        double predictedIntensity = (lastIntensity + psm.Query.spectrum.PrecursorIntensity) * 0.5;
                        intensityFactor = (psm.Query.spectrum.PrecursorIntensity - predictedIntensity) / predictedIntensity;// psm.Query.spectrum.PrecursorIntensity;
                    }
                    fragRatio.Add(psm, intensityFactor);
                }
            }
            return fragRatio;
        }

        public static List<string> GetFragments(Peptide peptide, int psmCharge, DBOptions dbOptions)
        {
            List<string> fragments = new List<string>();
            foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(peptide, psmCharge, dbOptions))
                {
                    if (fragment == match.fragment)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                    fragments.Add(fragment);
            }
            return fragments;
        }
        
        public static List<ProductMatch> GetCombinedSpectrum(DBOptions dbOptions, List<PeptideSpectrumMatch> psms, Peptide peptide, int psmCharge, Dictionary<double, int> DicOfCommonPM = null)
        {
            Dictionary<PeptideSpectrumMatch, double> DicOfPsmFactor = ComputeMsMsNormalizationFactors(psms);
            Dictionary<ProductMatch, double> DicOfProductMsMsFactor = new Dictionary<ProductMatch,double>();
            Dictionary<string, List<ProductMatch>> DicOfProducts = new Dictionary<string, List<ProductMatch>>();

            foreach (PeptideSpectrumMatch psm in psms)
            {
                foreach (ProductMatch match in psm.AllProductMatches)
                {
                    if (DicOfCommonPM == null || DicOfCommonPM.ContainsKey(match.theoMz))
                    {
                        string key = match.fragment + "|" + match.fragmentPos + "|" + match.charge;
                        if (!DicOfProducts.ContainsKey(key))
                            DicOfProducts.Add(key, new List<ProductMatch>());
                        DicOfProducts[key].Add(match);
                        DicOfProductMsMsFactor.Add(match, DicOfPsmFactor[psm]);
                    }
                }
            }

            double avgInt = 0;
            List<ProductMatch> products = new List<ProductMatch>();
            foreach (List<ProductMatch> matchList in DicOfProducts.Values)
            {
                ProductMatch newPM = new ProductMatch(matchList[0]);
                newPM.obsIntensity = 0;
                if (matchList.Count > 0)
                {
                    foreach (ProductMatch pm in matchList)
                        newPM.obsIntensity += pm.obsIntensity + pm.obsIntensity * DicOfProductMsMsFactor[pm];

                    newPM.obsIntensity /= (double)matchList.Count;
                }
                newPM.weight = matchList.Count * newPM.obsIntensity;
                avgInt += newPM.obsIntensity;
                products.Add(newPM);
            }
            avgInt /= (double) products.Count;

            //Add missed important fragments
            if (DicOfCommonPM != null)
            {
                foreach (double mz in DicOfCommonPM.Keys)
                {
                    bool found = false;
                    foreach (ProductMatch match in products)
                        if (match.theoMz == mz)
                            found = true;
                    if (!found)
                    {
                        ProductMatch newMatch = new ProductMatch();
                        newMatch.theoMz = mz;
                        newMatch.weight = 0;
                        newMatch.obsIntensity = 0;
                        foreach (PeptideSpectrumMatch psm in psms)
                        {
                            foreach (MsMsPeak peak in psm.Query.spectrum.Peaks)
                            {
                                if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(peak.MZ, mz, dbOptions.productMassTolerance.Units)) <= dbOptions.productMassTolerance.Value)
                                {
                                    newMatch.weight += 1;
                                    newMatch.obsIntensity += peak.Intensity + peak.Intensity * DicOfPsmFactor[psm];
                                }
                            }
                        }
                        if (newMatch.obsIntensity < avgInt * 0.05)
                            newMatch.obsIntensity = 0;
                        else
                            newMatch.obsIntensity /= (double)newMatch.weight;
                        newMatch.weight *= newMatch.obsIntensity;
                        products.Add(newMatch);
                    }
                }
            }
            return products;
        }
    }
}
