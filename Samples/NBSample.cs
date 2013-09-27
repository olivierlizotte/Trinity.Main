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
        public static long CreateFictiousSpectrum(List<double> frag2Ratio, int out2Ratio,
                                        List<double> frag3Ratio, int out3Ratio,
                                        List<double> frag4Ratio, int out4Ratio,
                                        List<double> frag5Ratio, int out5Ratio,
                                        List<long> fragments,
                                        ref long[] output)
        {
            long error = 0;
            for (int i = 0; i < frag2Ratio.Count; i++)
            {
                long sum = 0;
                sum += (long)(frag2Ratio[i] * out2Ratio);
                sum += (long)(frag3Ratio[i] * out3Ratio);
                sum += (long)(frag4Ratio[i] * out4Ratio);
                sum += (long)(frag5Ratio[i] * out5Ratio);
                output[i] = sum;
                error += Math.Abs(sum - fragments[i]);
            }
            return error;
        }

        public static long Overflow(List<List<double>> fragRatios, List<int> ratios,
                                        List<long> fragments)
        {
            long error = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                long sum = 0;
                for (int j = 0; j < ratios.Count; j++)
                    sum += (long)(fragRatios[j][i] * ratios[j]);
                if(sum > fragments[i])
                    error += sum - fragments[i];
            }
            return error;
        }

        public static long Underflow(List<List<double>> fragRatios, List<int> ratios,
                                        List<long> fragments)
        {
            long error = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                long sum = 0;
                for (int j = 0; j < ratios.Count; j++)
                    sum += (int)(fragRatios[j][i] * ratios[j]);
                if (sum < fragments[i])
                    error += fragments[i] - sum;
            }
            return error;
        }

        public static int FindLocalMaximumFlow(List<double> fragRatio, List<long> fragments, long sumOfIntensities)
        {
            int cumul = 0;
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

        public static long MaxFlow(List<List<double>> fragRatios,
                                    List<long> fragments,
                                ref List<List<double>> optimalSolutions)
        {
            //Lists must have same number of fragments, ordered in the same manner
            long sumOfIntensities = 0;
            for (int i = 0; i < fragments.Count; i++)
                sumOfIntensities += fragments[i];

            List<int> localFlows = new List<int>();
            List<int> localMaxFlows = new List<int>();
            foreach (List<double> fragmentRatio in fragRatios)
            {
                int maxFlow = FindLocalMaximumFlow(fragmentRatio, fragments, sumOfIntensities);
                localFlows.Add(maxFlow);
                localMaxFlows.Add(maxFlow);
            }

            long overError = Overflow(fragRatios, localFlows, fragments);
            long underError = Underflow(fragRatios, localFlows, fragments);
            int bestIndex = 0;

            int iterSize = 1;
            while (overError > 0 && iterSize < 10000)
            {
                bestIndex = -1;
                long smallestUnderError = long.MaxValue;
                for (int i = 0; i < fragRatios.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= iterSize;
                        long tmpErrorMinus = Underflow(fragRatios, localFlows, fragments);

                        if (tmpErrorMinus <= smallestUnderError && Overflow(fragRatios, localFlows, fragments) < overError)
                        {
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
            }
            Console.WriteLine("   Max Flow found!");
            optimalSolutions.Clear();

            List<double> newList = new List<double>();
            foreach (int localFlow in localFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            return Underflow(fragRatios, localFlows, fragments);
        }//*/

        public static long MaxFlowBruteForce(List<double> frag2Ratio,
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
            long smallestError = long.MaxValue;
            List<int[]> ratios = new List<int[]>();
            long[] outputs = new long[fragments.Count];
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
                            long error = CreateFictiousSpectrum(frag2Ratio, intensity2,
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

        public static Result Launch(bool restrain = false)
        {
            string outputDir = @"C:\_IRIC\DATA\Test\testNB\";
            string fastaFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\FEB13_2013\MSMS files\Peptide.fasta";//Yeast            

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

            //string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project__Group2_All.csv";//Group 2 (all)

            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project__Group3_All.csv";//Group 3 (all)
            
            //@"G:\Thibault\Olivier\MQ_vs_Morpheus\project.csv";//Yeast
            //@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JAN22_2013\_Project_FL_Single.csv";
            //G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUN27_2012\MR 4Rep DS\MassSense\_Test_ProjectFile_MF3.csv";
            //G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\MAR18_2013\ProjectFile_TestForProPheus.csv";
            Samples Project = new Samples(projectFile, 0);
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
            dbOptions.LoadSpectraIfFound = true;

            /*
            List<int> ace5Test = new List<int>(new int[] { 683,1456,0,1890,0,0,0,0,9005,0,2224,0,0,2224,0,4576,1643,1054,235,683 });
            List<int> ace8Test = new List<int>(new int[] { 532, 0, 0, 416, 1852, 2929, 0, 1286, 0, 0, 0, 0, 0, 2468, 0, 4539, 1855, 0, 409, 532 });
            List<int> ace12Test = new List<int>(new int[] { 575, 124, 0, 1724, 1824, 102, 110, 1264, 9028, 124, 2752, 2064, 396, 0, 2009, 55, 63, 0, 75, 575 });
            List<int> ace16Test = new List<int>(new int[] { 2108,6152,2110,2721,106,4819,891,468,125,6574,206,358,3131,71,574,248,679,287,565,2108 });
            List<int> mix = new List<int>();
            for (int i = 0; i < ace5Test.Count; i++)
                mix.Add((int)(ace5Test[i] * 0.1 + ace8Test[i] * 0.2 + ace12Test[i] * 0.6 + ace16Test[i] * 0.1));
            //Aimed     0.1,    0.2,    0.4,    0.3
            //At 50     0.1,    0.15,   0.37,   0.41

            //TestGroup2MaxFlow(dbOptions.OutputFolder + "MaxFlow_Group2.csv", mix);
            //*/

            Propheus propheus = new Propheus(dbOptions, Project);

            propheus.Preload();
            propheus.PrepareQueries();

            Result tmp = propheus.SearchLatestVersion(propheus.AllQueries, true);            
            tmp.WriteInfoToCsv(false);                      

            int charge = 2;
            Peptide peptide = null;   
            List<List<double>> items = new List<List<double>>();
            List<Sample> samples = new List<Sample>();
            foreach (Sample sample in Project)
            {
                List<PeptideSpectrumMatch> psmList = new List<PeptideSpectrumMatch>();
                foreach (Precursor precursor in tmp.matchedPrecursors)
                    if (precursor.Charge == charge && precursor.sample == sample)
                    {
                        foreach (PeptideSpectrumMatch psm in precursor.psms)
                            if (sample.nameColumn.CompareTo(psm.Peptide.Sequence) == 0)
                            {
                                peptide = psm.Peptide;
                                psmList.Add(psm);
                            }
                    }
                List<double> tmpList = tmp.GetFragments(psmList, peptide, charge, false, false);                
                items.Add(tmpList);
                samples.Add(sample);
            }

            double FragmentStepSize = 1000;
            double RatiosStepSize = 1000;
            //double 
            foreach (List<double> list in items)
                NormalizeList(list, RatiosStepSize);//*/

            //Build InSilico Mixed Spectrum
            List<long> mix = new List<long>();
            for (int i = 0; i < items[0].Count; i++)
            {
                double cumul = 0;
                double iter = 5;
                foreach (List<double> list in items)
                {
                    cumul += list[i] * iter;
                    iter *= 2;
                }
                mix.Add((long)(cumul * FragmentStepSize));
            }

            List<List<double>> optimalSolutions = new List<List<double>>();
            long error = MaxFlow(items, mix, ref optimalSolutions);            
            //int error = MaxFlowBruteForce(items[0], items[1], items[2], items[3], mix, ref optimalSolutions);            
            foreach (List<double> solution in optimalSolutions)
            {
                Console.WriteLine("Max Flow computed (error of " + error / (FragmentStepSize * RatiosStepSize) + ")");
                for(int i = 0; i < samples.Count; i++)
                    Console.WriteLine("     " + samples[i].nameColumn + " -> " + solution[i] / FragmentStepSize + "");
            }
            Console.WriteLine("Number of solutions : " + optimalSolutions.Count);
            return tmp;
        }
    }
}
