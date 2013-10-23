﻿/*
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
    public class IsomericPeptideQuantifier
    {
        public static double ComputeOverflow(List<List<double>> fragRatios, List<int> ratios, List<long> fragments)
        {
            double error = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                double sum = 0;
                for (int j = 0; j < ratios.Count; j++)
                    sum += fragRatios[j][i] * ratios[j];
                if(sum > fragments[i])
                    error += sum - fragments[i];
            }
            return error;
        }

        public static double ComputeUnderflow(List<List<double>> fragRatios, List<int> ratios,
                                        List<long> fragments)
        {
            double error = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                double sum = 0;
                for (int j = 0; j < ratios.Count; j++)
                    sum += fragRatios[j][i] * ratios[j];
                if (sum < fragments[i])
                    error += fragments[i] - sum;
            }
            return error;
        }

        public static int FindLocalMaximumFlow(List<double> fragRatio, List<long> fragments, long sumOfIntensities)
        {
            int cumul = 1;
            while (cumul <= sumOfIntensities)
            {
                for (int i = 0; i < fragments.Count; i++)
                {
                    if (fragRatio[i] * cumul > fragments[i])
                        return cumul;
                }
                cumul++;
            }
            return cumul;
        }

        public static double MaxFlow(List<List<double>> fragRatios,
                                    List<long> fragments,
                                ref List<List<double>> optimalSolutions)
        {
            //Lists must have same number of fragments, ordered in the same manner
            long sumOfIntensities = 0;
            for (int i = 0; i < fragments.Count; i++)
                sumOfIntensities += fragments[i];

            List<int> localFlows = new List<int>();
            foreach (List<double> fragmentRatio in fragRatios)
                localFlows.Add(FindLocalMaximumFlow(fragmentRatio, fragments, sumOfIntensities));

            double overError    = ComputeOverflow(fragRatios, localFlows, fragments);
            double underError   = ComputeUnderflow(fragRatios, localFlows, fragments);
            int bestIndex = 0;

            int iterSize = 1;
            double bestOverallError = double.MaxValue;
            List<int> bestLocalFlows = new List<int>();

            while (overError > 0 && iterSize < 10000)
            {
                bestIndex = -1;
                double smallestUnderError = double.MaxValue;
                double smallestOverError = overError;
                for (int i = 0; i < fragRatios.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= iterSize;
                        double tmpErrorMinus = ComputeUnderflow(fragRatios, localFlows, fragments);

                        if (tmpErrorMinus < smallestUnderError// && Overflow(fragRatios, localFlows, fragments) < overError)
                            || ( tmpErrorMinus == smallestUnderError && ComputeOverflow(fragRatios, localFlows, fragments) < smallestOverError))
                        {
                            smallestOverError  = ComputeOverflow(fragRatios, localFlows, fragments);
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

                overError = ComputeOverflow(fragRatios, localFlows, fragments);
                underError = ComputeUnderflow(fragRatios, localFlows, fragments);
                if (overError + underError < bestOverallError)
                {
                    bestLocalFlows = new List<int>(localFlows);
                    bestOverallError = overError + underError;
                }
            }
            optimalSolutions.Clear();

            List<double> newList = new List<double>();
            foreach (int localFlow in localFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            newList = new List<double>();
            foreach (int localFlow in bestLocalFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            return ComputeUnderflow(fragRatios, localFlows, fragments);
        }//*/

        public static DBOptions GetDBOptions(bool loadFromRaw)
        {
            string outputDir = @"C:\_IRIC\DATA\Test\testNB\";
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
            dbOptions.fragments.Add(new FragmentA());
            dbOptions.fragments.Add(new FragmentB());
            dbOptions.fragments.Add(new FragmentC());
            dbOptions.fragments.Add(new FragmentX());
            dbOptions.fragments.Add(new FragmentY());
            dbOptions.fragments.Add(new FragmentZ());

            //ClusterOptions clusterOptions = new ClusterOptions(Project, outputDir, 5, true, 90, true);//TODO validate its in seconds for all file types

            dbOptions.SaveMS1Peaks = true;
            dbOptions.SaveMSMSPeaks = true;
            dbOptions.LoadSpectraIfFound = !loadFromRaw;

            dbOptions.NbPSMToKeep = 64;
            return dbOptions;
        }
        
        public static void MaxFlowThis()//string projectSingleInjections, string projectMixed)
        {
            Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Varied_19Oct.csv", 0);

            string baseSeq = "";
            List<List<ProductMatch>> ratios = new List<List<ProductMatch>>();
            List<double> TrackIntensity = new List<double>();
            List<double> NormalizeFactor = new List<double>();
            CreateVirtualSpectrumFromSpikedPeptides(ProjectRatios, ref baseSeq, false, false, 0.01, ref ratios, ref TrackIntensity, ref NormalizeFactor);

            List<string> ratioNames = new List<string>();
            foreach (Sample sample in ProjectRatios)
                ratioNames.Add(sample.nameColumn);

            //Get PSMs
            DBOptions dbOptions = GetDBOptions(false);
            Propheus propheus = new Propheus(dbOptions, ProjectMixed);

            propheus.Preload(false, false);
            propheus.PrepareQueries();

            Result tmp = propheus.SearchLatestVersion(propheus.AllQueries, false);

            int iterPsm = 0;
            double bestScore = double.MaxValue;
            vsCSVWriter writerCumul = new vsCSVWriter(dbOptions.OutputFolder + "CumulRatios.csv");
            foreach(Sample sample in ProjectMixed)
            {
                string lineCumulRatio = sample.sSDF;
                //PeptideSpectrumMatch bestPsm = null;
                Console.WriteLine("Sample " + sample.sSDF);
                /*foreach (Query query in tmp.queries)
                {
                    if (query.sample == sample && query.spectrum.ScanNumber == 4635)
                    {
                        foreach (PeptideSpectrumMatch psm in query.psms)
                        {
                            bool keep = false;
                            foreach (string name in ratioNames)
                                if (name.CompareTo(psm.Peptide.Sequence) == 0)
                                    keep = true;
                            if (keep && psm.Peptide.BaseSequence.CompareTo(baseSeq) == 0 && query.precursor.Charge == 2 && psm.MatchingProducts > 4)
                            {
                                if (bestPsm == null || 
                                    psm.MatchingIntensityFraction > bestPsm.MatchingIntensityFraction)
                                    //query.spectrum.PrecursorIntensity > bestPsm.Query.spectrum.PrecursorIntensity)// ||
                                    //(query.spectrum.PrecursorIntensity > bestPsm.Query.spectrum.PrecursorIntensity && psm.MatchingIntensityFraction > bestPsm.MatchingIntensityFraction))
                                    bestPsm = psm;
                            }
                        }
                    }
                }
                if(bestPsm != null)
                {//*/
                Dictionary<ProductSpectrum, bool> doneSpectrum = new Dictionary<ProductSpectrum, bool>();
                Dictionary<int, string> dicOfResults = new Dictionary<int, string>();
                List<double> sumOfRatio = new List<double>();
                foreach (double intensity in TrackIntensity)
                    sumOfRatio.Add(0);

                foreach (Query query in tmp.queries)
                {
                    if (query.sample == sample)
                    {
                        foreach (PeptideSpectrumMatch psm in query.psms)
                        {
                            bool keep = false;
                            foreach (string name in ratioNames)
                                if (name.CompareTo(psm.Peptide.Sequence) == 0)
                                    keep = true;
                            
                            if (keep && !doneSpectrum.ContainsKey(query.spectrum) && psm.Peptide.BaseSequence.CompareTo(baseSeq) == 0 && query.precursor.Charge == 2 && psm.MatchingProducts > 4)
                            {
                                doneSpectrum.Add(query.spectrum, true);
                                iterPsm++;

                                //List<double> capacity = new List<double>();
                                //foreach (double fragment in tmp.GetFragments(psm, psm.Peptide, 2))
                                //    capacity.Add(fragment);

                                double overFlow = 0;
                                double underFlow = 0;
                                double percentError = 0;
                                List<double> finalRatios = MaxFlowFromSpectrum(ratios, ratioNames, 1000000, query.spectrum.Peaks, dbOptions.productMassTolerance, ref overFlow, ref underFlow, ref percentError);
                                bool AllThere = true;
                                foreach (double dbl in finalRatios)
                                    if (dbl == 0)
                                        AllThere = false;

                                double normSumRatio = 0;
                                for (int i = 0; i < finalRatios.Count; i++)
                                    normSumRatio += finalRatios[i] / NormalizeFactor[i];

                                double sumRatio = 0;
                                for (int i = 0; i < finalRatios.Count; i++)
                                {
                                    sumRatio += finalRatios[i];
                                }

                                if (AllThere || percentError < bestScore)//percentError < 50)
                                {
                                    if (percentError < bestScore)
                                        bestScore = percentError;
                                    Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");
                                    Console.WriteLine("Max Flow computed (error of " + (float)percentError + "%)");
                                    for (int i = 0; i < ratioNames.Count; i++)
                                        Console.WriteLine(ratioNames[i] + " -> [ " + ((float)(NormalizeFactor[i] * finalRatios[i] / normSumRatio)) * 100 + " ; " + ((float)(finalRatios[i] / sumRatio)) * 100 + " ; " + (float)finalRatios[i] + " ; " + (float)(NormalizeFactor[i] * finalRatios[i]) + " ; " + (finalRatios[i] / sumRatio) * query.precursor.Track.INTENSITY * (500 / TrackIntensity[i]) + " ]");
                                    Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");
                                }
                                string strRatios = query.spectrum.RetentionTimeInMin.ToString() + "," + psm.Query.spectrum.PrecursorIntensity;
                                //for (int i = 0; i < finalRatios.Count; i++)
                                //    strRatios += "," + finalRatios[i];
                                //for (int i = 0; i < finalRatios.Count; i++)
                                //    strRatios += "," + finalRatios[i] / NormalizeFactor[i];
                                for (int i = 0; i < finalRatios.Count; i++)
                                {
                                    sumOfRatio[i] += (finalRatios[i] * NormalizeFactor[i]);
                                    strRatios += "," + finalRatios[i] * NormalizeFactor[i];
                                }
                                for (int i = 0; i < finalRatios.Count; i++)
                                    strRatios += "," + (finalRatios[i] * NormalizeFactor[i]) * query.spectrum.TotalIntensity;
                                for (int i = 0; i < finalRatios.Count; i++)
                                {
                                    //sumOfRatio[i] +=   (finalRatios[i] * NormalizeFactor[i]) * query.precursor.Track.INTENSITY * (500 / TrackIntensity[i]);
                                    strRatios += "," + (finalRatios[i] * NormalizeFactor[i]) * query.precursor.Track.INTENSITY * (500 / TrackIntensity[i]);
                                }
                                for (int i = 0; i < finalRatios.Count; i++)
                                    strRatios += "," + ((finalRatios[i] * NormalizeFactor[i]) * query.spectrum.TotalIntensity) * query.precursor.Track.INTENSITY * (500 / TrackIntensity[i]);
                                dicOfResults.Add(query.spectrum.ScanNumber, strRatios);
                            }
                        }
                    }
                }
                vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + sample.nameColumn + "_Ratios.csv");
                List<int> scans = new List<int>(dicOfResults.Keys);
                scans.Sort();
                foreach (int scan in scans)
                    writer.AddLine(dicOfResults[scan]);
                writer.writeToFile();
                foreach (double n in sumOfRatio)
                    lineCumulRatio += "," + n;
                writerCumul.AddLine(lineCumulRatio);
            }
            writerCumul.writeToFile();
        }

        private static double ComputeMaxFlow(List<List<ProductMatch>> spikedMatches,
                                    List<MsMsPeak> mixedSpectrum, MassTolerance tolerance,
                                ref List<List<double>> optimalSolutions,
                                ref double percentError,
                                ref List<long> average)
        {
            //Lists must have same number of fragments, ordered in the same manner
            double sumOfIntensities = 0;
            for (int i = 0; i < mixedSpectrum.Count; i++)
                sumOfIntensities += mixedSpectrum[i].Intensity;

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

            while (overError > 0 && iterSize < 10000)
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

            virtualSpectrum = BuildVirtualSpectrum(spikedMatches, bestLocalFlows, mixedFragDic);
            //Compute expected error in percentage
            percentError = (ComputeOverflow(virtualSpectrum, mixedFragDic) / sumOfIntensities) * 100.0;
            
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

        public static List<double> MaxFlowFromSpectrum(List<List<ProductMatch>> ratiosToFit, List<string> ratioNames, int precision, GraphML_List<MsMsPeak> capacity, MassTolerance tolerance, 
            ref double overFlow, ref double underFlow, ref double errorInPercent)
        {
            //double 
            /*
            foreach (List<double> list in ratiosToFit)
                NormalizeList(list, precision);//*/

            List<List<double>> solutions = new List<List<double>>();
            List<long> average = new List<long>();
            List<MsMsPeak> expandedCapacity = new List<MsMsPeak>();
            foreach (MsMsPeak peak in capacity)
                expandedCapacity.Add(new MsMsPeak(peak.MZ, peak.Intensity * precision, peak.Charge));

            double error = ComputeMaxFlow(ratiosToFit, expandedCapacity, tolerance, ref solutions, ref errorInPercent, ref average);

            double sumOfIntensities = 0;
            foreach (MsMsPeak peak in expandedCapacity)
                sumOfIntensities += peak.Intensity;
            /*
            Console.WriteLine(" -=+ Error cumulated : " + (float)(errorCumul / (double)sumOfIntensities) * 100 + "% +=- ");
            
            foreach (List<double> solution in solutions)
            {
                Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");
                Console.WriteLine("Max Flow computed (error of " + error / (double)(precision * precision) + ")");
                for (int i = 0; i < ratioNames.Count; i++)
                    Console.WriteLine("     " + ratioNames[i] + " -> " + solution[i] / (double) precision + "");
            }
            Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");
            Console.WriteLine("Number of solutions : " + solutions.Count);
            //*/
            overFlow = 0;
            underFlow = error;

            List<double> result = new List<double>();            
            //foreach(double val in solutions[1])
            //    result.Add(val / (double)precision);
            foreach (double val in average)
                result.Add(val / (double)precision);
            return result;
        }

        public static void CreateVirtualSpectrumFromSpikedPeptides(Samples Project, ref string BaseSequence, bool loadFromRaw, bool bestPSMOnly, double minRatioIntensity,
                    ref List<List<ProductMatch>> ratios, ref List<double> TrackIntensity, ref List<double> Normalizor)
        {
            DBOptions dbOptions = GetDBOptions(loadFromRaw);
            Propheus propheus = new Propheus(dbOptions, Project);

            propheus.Preload(false, false);
            propheus.PrepareQueries();

            Result tmp = propheus.SearchLatestVersion(propheus.AllQueries, false);
            tmp.WriteInfoToCsv(false);

            int charge = 2;
            Peptide peptide = null;
            ratios.Clear();

            double averagePrecursorIntensity = 0;

            foreach (Sample sample in Project)
            {
                double precursorMaxIntensity = 0;
                List<PeptideSpectrumMatch> psmList = new List<PeptideSpectrumMatch>();
                PeptideSpectrumMatch bestPsm = null;
                foreach (Query query in tmp.queries)
                {
                    if (query.precursor.Charge == charge && query.precursor.sample == sample)
                    {
                        foreach (PeptideSpectrumMatch psm in query.psms)
                            if (sample.nameColumn.CompareTo(psm.Peptide.Sequence) == 0 && psm.MatchingProducts > 4)//.ProbabilityScore() > 0.020)// && psm.ProbabilityScore() > bestPsmScore)
                            {
                                peptide = psm.Peptide;
                                BaseSequence = peptide.BaseSequence;
                                psmList.Add(psm);
                                //if (!query.precursor.Track.Invented)
                                //{
                                //    nbPrec++;
                                if (bestPsm == null || psm.MatchingIntensity < bestPsm.MatchingIntensity)
                                {
                                    bestPsm = psm;
                                    if(bestPSMOnly)
                                        precursorMaxIntensity = query.spectrum.TotalIntensity;
                                }
                                else if (!bestPSMOnly && precursorMaxIntensity < query.precursor.Track.INTENSITY)
                                    {
                                        precursorMaxIntensity = query.precursor.Track.INTENSITY;
                                    }
                                
                                //}
                                //else
                                //    nbPrec = nbPrec;
                            }
                    }
                }
                if (bestPSMOnly)
                {
                    psmList.Clear();
                    psmList.Add(bestPsm);
                }
                else
                {
                    psmList.Sort(PeptideSpectrumMatches.CompareMatchingIntensity);
                    if (psmList.Count > 20)
                        psmList.RemoveRange(20, psmList.Count - 20);
                }
                List<ProductMatch> tmpList = tmp.GetCommonSpectrum(psmList, peptide, charge);
                //double factor = 1 / (double)psmList.Count;
                /*
                double HighestIntentity= 0;
                foreach (ProductMatch elem in tmpList)
                    if (elem.obsIntensity > HighestIntentity)
                        HighestIntentity = elem.obsIntensity;

                factor = 10 * (1.0 / HighestIntentity);
                for(int i = 0; i < tmpList.Count; i++)
                    tmpList[i].obsIntensity *= factor;//*/
                //Normalizor.Add(factor);
                averagePrecursorIntensity += precursorMaxIntensity;

                ratios.Add(tmpList);
                TrackIntensity.Add(precursorMaxIntensity);// / (double)nbPrec);
                tmp.ExportFragmentIntensitiesForAllPSM(psmList, peptide, charge, dbOptions.OutputFolder + sample.nameColumn + ".csv");
            }
            averagePrecursorIntensity /= (double)Project.Count;
            for (int i = 0; i < Project.Count; i++)
            {
                double factor = averagePrecursorIntensity / TrackIntensity[i];
                for (int j = 0; j < ratios[i].Count; j++)
                    ratios[i][j].obsIntensity *= factor;//*/
                Normalizor.Add(factor);
            }
        }
    }
}
