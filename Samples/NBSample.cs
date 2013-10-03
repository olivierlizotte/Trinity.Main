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
    public class NBSample
    {
        public static double CreateFictiousSpectrum(List<double> frag2Ratio, int out2Ratio,
                                        List<double> frag3Ratio, int out3Ratio,
                                        List<double> frag4Ratio, int out4Ratio,
                                        List<double> frag5Ratio, int out5Ratio,
                                        List<long> fragments,
                                        ref double[] output)
        {
            double error = 0;
            for (int i = 0; i < frag2Ratio.Count; i++)
            {
                double sum = 0;
                sum += frag2Ratio[i] * out2Ratio;
                sum += frag3Ratio[i] * out3Ratio;
                sum += frag4Ratio[i] * out4Ratio;
                sum += frag5Ratio[i] * out5Ratio;
                output[i] = sum;
                error += Math.Abs(sum - fragments[i]);
            }
            return error;
        }

        public static double Overflow(List<List<double>> fragRatios, List<int> ratios, List<long> fragments)
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

        public static double Underflow(List<List<double>> fragRatios, List<int> ratios,
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

            double overError = Overflow(fragRatios, localFlows, fragments);
            double underError = Underflow(fragRatios, localFlows, fragments);
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
                        double tmpErrorMinus = Underflow(fragRatios, localFlows, fragments);

                        if (tmpErrorMinus < smallestUnderError// && Overflow(fragRatios, localFlows, fragments) < overError)
                            || ( tmpErrorMinus == smallestUnderError && Overflow(fragRatios, localFlows, fragments) < smallestOverError))
                        {
                            smallestOverError  = Overflow(fragRatios, localFlows, fragments);
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

                overError = Overflow(fragRatios, localFlows, fragments);
                underError = Underflow(fragRatios, localFlows, fragments);
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

            return Underflow(fragRatios, localFlows, fragments);
        }//*/

        public static DBOptions GetDBOptions(bool loadFromRaw)
        {
            string outputDir = @"C:\_IRIC\DATA\Test\testNB\";
            string fastaFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\FEB13_2013\MSMS files\Peptide.fasta";
            DBOptions dbOptions = new DBOptions(fastaFile);
            dbOptions.precursorMassTolerance = new MassTolerance(8/*8*//*8withoutisotopes*/, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(0.5/*0.034*//*without isotopes*/, MassToleranceUnits.Da);//0.034 is a 60 000 resolution over 2000 range in mz
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

            dbOptions.addFragmentLoss = false;
            dbOptions.addFragmentMods = true;
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

            return dbOptions;
        }

        public static double MaxFlowBruteForce(List<double> frag2Ratio,
                                List<double> frag3Ratio,
                                List<double> frag4Ratio,
                                List<double> frag5Ratio,
                                List<long> fragments,
                                ref List<List<double>> optimalSolutions)
        {
            //Lists must have same number of fragments, ordered in the same manner
            int stepSize = 10;// 100;
            long sum = 0;
            long sumOfIntensities = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                sumOfIntensities += fragments[i];
                fragments[i] = fragments[i] / stepSize;
                sum += fragments[i];
            }
            sum += stepSize;
            double smallestError = double.MaxValue;
            List<int[]> ratios = new List<int[]>();
            double[] outputs = new double[fragments.Count];
            int localMax2 = FindLocalMaximumFlow(frag2Ratio, fragments, sum);
            int localMax3 = FindLocalMaximumFlow(frag3Ratio, fragments, sum);
            int localMax4 = FindLocalMaximumFlow(frag4Ratio, fragments, sum);
            int localMax5 = FindLocalMaximumFlow(frag5Ratio, fragments, sum);
            for (int intensity2 = 0; intensity2 < localMax2; intensity2++)
            {
                double max3 = sum - intensity2;
                if (localMax3 < max3)
                    max3 = localMax3;
                for (int intensity3 = 0; intensity3 < max3; intensity3++)
                {
                    double max4 = sum - intensity2 - intensity3;
                    if (localMax4 < max4)
                        max4 = localMax4;
                    for (int intensity4 = 0; intensity4 < max4; intensity4++)
                    {
                        double max5 = sum - intensity2 - intensity3 - intensity4;
                        if (localMax5 < max5)
                            max5 = localMax5;
                        for (int intensity5 = 0; intensity5 < max5; intensity5++)
                        {
                            double error = CreateFictiousSpectrum(frag2Ratio, intensity2,
                                                    frag3Ratio, intensity3,
                                                    frag4Ratio, intensity4,
                                                    frag5Ratio, intensity5,
                                                    fragments,
                                                    ref outputs);
                            
                            if (error <= smallestError)
                            {
                                bool isOk = true;
                                if (error == 0)
                                    isOk = true;
                                if (isOk)
                                {
                                    if (error != smallestError)
                                        ratios = new List<int[]>();
                                    ratios.Add(new int[4] { intensity2, intensity3, intensity4, intensity5 });
                                    smallestError = error;
                                }
                            }
                        }
                    }
                }

                Console.Write("\r{0}%   ", ((100 * intensity2) / localMax2));
            }
            Console.WriteLine("   Max Flow found!");
            optimalSolutions.Clear();
            for (int i = 0; i < ratios.Count; i++)
            {
                List<double> newList = new List<double>();
                for (int j = 0; j < ratios[i].Length; j++)
                    newList.Add(stepSize * (ratios[i][j] / (double)sumOfIntensities));
                optimalSolutions.Add(newList);
            }

            return smallestError * stepSize;// + sum;
        }

        public static List<double> BuildFictitiousRatios(int nbItems, Random rd)
        {
            List<double> list = new List<double>();
            double sum = 0;
            for (int i = 0; i < nbItems - 1; i++)
            {
                double ratio = (1 - sum) * rd.NextDouble();
                list.Add(ratio);
                sum += ratio;
            }
            list.Add(1 - sum);
            return list;
        }

        public static void NormalizeList(List<double> list, double dividend)
        {
            double sum = 0.0;
            foreach (double item in list)
                sum += item;
            sum /= dividend;
            for (int i = 0; i < list.Count; i++)
                list[i] /= sum;
        }

        public static void MaxFlowThis()//string projectSingleInjections, string projectMixed)
        {
            Samples ProjectRatios = new Samples(@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project__Group2_All.csv", 0);//Group 2 (all)            
            Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_OneRaw.csv", 0);

            string baseSeq = "";
            List<List<double>> ratios = BuildPeptideRatios(ProjectRatios, ref baseSeq, false, true, 0.01);
            
            double RatiosStepSize = 1000000;
            foreach (List<double> list in ratios)
                NormalizeList(list, RatiosStepSize);//*/
            
            //Get PSMs
            DBOptions dbOptions = GetDBOptions(false);
            dbOptions.productMassTolerance.Value = 0.05;
            Propheus propheus = new Propheus(dbOptions, ProjectMixed);

            propheus.Preload();
            propheus.PrepareQueries();

            Result tmp = propheus.SearchLatestVersion(propheus.AllQueries, true);

            List<PeptideSpectrumMatch> tmpList = new List<PeptideSpectrumMatch>();
            foreach (Query query in tmp.queries)
            {
                foreach (PeptideSpectrumMatch psm in query.psms)
                {
                    if (psm.Peptide.BaseSequence.CompareTo(baseSeq) == 0 && query.precursor.Charge == 2)
                    {
                        tmpList.Add(psm);
                        break;
                    }
                }
            }
            
            List<long> mix = new List<long>();
            double FragmentStepSize = 1000;
            foreach (double fragment in tmp.GetFragments(tmpList, tmpList[0].Peptide, 2, false, false, 0.0001))
                mix.Add((long)(fragment * FragmentStepSize * RatiosStepSize));

            List<List<double>> optimalSolutions = new List<List<double>>();
            double error = MaxFlow(ratios, mix, ref optimalSolutions);

            bool isEmpty = true;
            foreach (List<double> listOp in optimalSolutions)
                foreach (double dbl in listOp)
                    if (dbl > 0)
                        isEmpty = false;
            if (!isEmpty)
            {
                //Compute average
                List<double> average = new List<double>();
                for (int i = 0; i < optimalSolutions[0].Count; i++)
                {
                    double sum = 0.0;
                    foreach (List<double> solution in optimalSolutions)
                        sum += solution[i];
                    double avg = sum / (double)optimalSolutions.Count;
                    average.Add(avg);
                }
                optimalSolutions.Add(average);

                //Compute expected error in percentage
                double errorCumul = 0.0;

                for (int k = 0; k < ratios[0].Count; k++)
                {
                    double peakSum = 0;
                    for (int i = 0; i < ratios.Count; i++)
                        peakSum += ratios[i][k] * optimalSolutions[0][i];
                    errorCumul += Math.Abs(mix[k] - peakSum);
                }
                double sumOfIntensities = 0;
                foreach (long item in mix)
                    sumOfIntensities += item;
                //errorCumul /= (double)average.Count;
                Console.WriteLine(" -=+ Error cumulated : " + (float)(errorCumul / (double)sumOfIntensities) * 100 + " +=- ");
                List<Sample> samples = new List<Sample>();
                foreach (Sample sample in ProjectRatios)
                    samples.Add(sample);
                //int error = MaxFlowBruteForce(items[0], items[1], items[2], items[3], mix, ref optimalSolutions);            
                foreach (List<double> solution in optimalSolutions)
                {
                    Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");
                    Console.WriteLine("Max Flow computed (error of " + error / (FragmentStepSize * RatiosStepSize) + ")");
                    for (int i = 0; i < samples.Count; i++)
                        Console.WriteLine("     " + samples[i].nameColumn + " -> " + solution[i] / FragmentStepSize + "");
                }
                Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");

                Console.WriteLine("Number of solutions : " + optimalSolutions.Count);            
            }
        }

        public static List<List<double>> BuildPeptideRatios(Samples Project, ref string BaseSequence, bool loadFromRaw, bool bestPSMOnly, double minRatioIntensity)
        {
            DBOptions dbOptions = GetDBOptions(loadFromRaw);
            Propheus propheus = new Propheus(dbOptions, Project);

            propheus.Preload();
            propheus.PrepareQueries();

            Result tmp = propheus.SearchLatestVersion(propheus.AllQueries, true);
            tmp.WriteInfoToCsv(false);

            int charge = 2;
            Peptide peptide = null;
            List<List<double>> items = new List<List<double>>();
            foreach (Sample sample in Project)
            {
                List<PeptideSpectrumMatch> psmList = new List<PeptideSpectrumMatch>();
                double bestPsmScore = 0.0;
                PeptideSpectrumMatch bestPSM = null;
                foreach (Query query in tmp.queries)
                    if (query.precursor.Charge == charge && query.precursor.sample == sample)
                    {
                        foreach (PeptideSpectrumMatch psm in query.psms)
                            if (sample.nameColumn.CompareTo(psm.Peptide.Sequence) == 0 && psm.ProbabilityScore() > 0)// && psm.ProbabilityScore() > bestPsmScore)
                            {
                                peptide = psm.Peptide;
                                BaseSequence = peptide.BaseSequence;
                                if (psm.Query.precursor.Track.INTENSITY/*.ProbabilityScore()*/ > bestPsmScore)
                                {
                                    bestPsmScore = psm.Query.precursor.Track.INTENSITY;
                                    bestPSM = psm;
                                }
                                psmList.Add(psm);
                            }
                    }

                if (bestPSMOnly && bestPSM != null)
                {
                    //TODO Check witch one is better : Average, Sum or Best ProbabilityScore
                    psmList.Clear();
                    psmList.Add(bestPSM);
                }
                List<double> tmpList = tmp.GetFragments(psmList, peptide, charge, false, false, minRatioIntensity);
                items.Add(tmpList);

                tmp.ExportFragmentIntensitiesForAllPSM(psmList, peptide, charge, dbOptions.OutputFolder + sample.nameColumn + ".csv");
            }
            return items;
        }
        
        public static void Launch(bool restrain = false)
        {
            /*
            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project_1.csv";//Group 1
            string strToAim = "GK(propionylation of K)GGK(propionylation of K)GLGK(propionylation of K)GGAK(propionylation of K)R";
            double mzAIM = 747.94000244140625;//*/
            /*
            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project_2.csv";//Group 2
            string strToAim = "GK(acetylation of K)GGK(propionylation of K)GLGK(propionylation of K)GGAK(propionylation of K)R";
            double mzAIM = 747.93000244140625;//*/
            /*
            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project_3.csv";//Group 2
            string strToAim = "GK(propionylation of K)GGK(acetylation of K)GLGK(propionylation of K)GGAK(propionylation of K)R";
            double mzAIM = 747.93000244140625;//*/
            /*
            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project_4.csv";//Group 2
            string strToAim = "GK(propionylation of K)GGK(propionylation of K)GLGK(acetylation of K)GGAK(propionylation of K)R";
            double mzAIM = 747.93000244140625;//*/
            /*
            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project_5.csv";//Group 2
            string strToAim = "GK(propionylation of K)GGK(propionylation of K)GLGK(propionylation of K)GGAK(acetylation of K)R";
            double mzAIM = 747.93000244140625;//*/

            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project__Group2_All.csv";//Group 2 (all)

            //string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project__Group3_All.csv";//Group 3 (all)
            Samples Project = new Samples(projectFile, 0);

            string baseSeq = "";
            List<List<double>> items = BuildPeptideRatios(Project, ref baseSeq, false, false, 0.01);

            double FragmentStepSize = 1000;
            double RatiosStepSize = 1000;
            //double 
            foreach (List<double> list in items)
                NormalizeList(list, RatiosStepSize);//*/

            //Build InSilico Mixed Spectrum
            Random rd = new Random();
            List<long> mix = new List<long>();
            for (int i = 0; i < items[0].Count; i++)
            {
                double cumul = 0;
                //double iter = 640;
                double iter = 5;
                foreach (List<double> list in items)
                {
                    cumul += list[i] * iter;
                    //iter /= 2.0;
                    iter *= 2.0;
                }
                double delta = 1 + (rd.NextDouble() / 10.0) - 0.05;// * FragmentStepSize * RatiosStepSize);
                mix.Add((long)(cumul * FragmentStepSize * delta));//Pollute system
            }

            List<List<double>> optimalSolutions = new List<List<double>>();
            double error = MaxFlow(items, mix, ref optimalSolutions);   

            //Compute average
            List<double> average = new List<double>();
            for(int i = 0; i < optimalSolutions[0].Count; i++)
            {
                double sum = 0.0;
                foreach (List<double> solution in optimalSolutions)                
                    sum += solution[i];
                double avg = sum / (double) optimalSolutions.Count;
                average.Add(avg);
            }
            optimalSolutions.Add(average);

            //Compute expected error in percentage
            double errorCumul = 0.0;
            
            for(int k = 0; k < items[0].Count; k++)
            {
                double peakSum = 0;
                for(int i = 0; i < items.Count; i++)
                    peakSum += items[i][k] * optimalSolutions[0][i];
                errorCumul += Math.Abs(mix[k] - peakSum);

                //for(int k = 0; k < mix.Count; k++)
                //    errorPeak += Math.Abs(optimalSolutions[0][i] * mix[k] - optimalSolutions[1][i] * mix[k]);
                //foreach (List<double> solution in optimalSolutions)
                //    if(average[i] > 0)
                //        errorPeak += Math.Abs(solution[i] * mix) / average[i];
                //errorCumul += errorPeak;
            }
            double sumOfIntensities = 0;
            foreach (long item in mix)
                sumOfIntensities += item;
            //errorCumul /= (double)average.Count;
            Console.WriteLine(" -=+ Error cumulated : " + (float)(errorCumul / (double) sumOfIntensities) * 100 + " +=- ");
            List<Sample> samples = new List<Sample>();
            foreach (Sample sample in Project)
                samples.Add(sample);
            //int error = MaxFlowBruteForce(items[0], items[1], items[2], items[3], mix, ref optimalSolutions);            
            foreach (List<double> solution in optimalSolutions)
            {   
                Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");
                Console.WriteLine("Max Flow computed (error of " + error / (FragmentStepSize * RatiosStepSize) + ")");
                for(int i = 0; i < samples.Count; i++)
                    Console.WriteLine("     " + samples[i].nameColumn + " -> " + solution[i] / FragmentStepSize + "");
            }
            Console.WriteLine(" ----------------- --------------- ----------------- ------------------ ---------------- ");

            Console.WriteLine("Number of solutions : " + optimalSolutions.Count);
            //return tmp;
        }
    }
}
