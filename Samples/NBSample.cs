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
        public static int CreateFictiousSpectrum(  List<double> frag2Ratio, int out2Ratio,
                                        List<double> frag3Ratio, int out3Ratio,
                                        List<double> frag4Ratio, int out4Ratio,
                                        List<double> frag5Ratio, int out5Ratio,
                                        List<int> fragments,
                                        ref int[] output)
        {
            int error = 0;
            for (int i = 0; i < frag2Ratio.Count; i++)
            {
                int sum = 0;
                sum += (int)(frag2Ratio[i] * out2Ratio);
                sum += (int)(frag3Ratio[i] * out3Ratio);
                sum += (int)(frag4Ratio[i] * out4Ratio);
                sum += (int)(frag5Ratio[i] * out5Ratio);
                output[i] = sum;
                error += Math.Abs(sum - fragments[i]);
            }
            return error;
        }

        public static int Overflow(List<List<double>> fragRatios, List<int> ratios,
                                        List<int> fragments)
        {
            int error = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                int sum = 0;
                for (int j = 0; j < ratios.Count; j++)
                    sum += (int)(fragRatios[j][i] * ratios[j]);
                if(sum > fragments[i])
                    error += sum - fragments[i];
            }
            return error;
        }

        public static int Underflow(List<List<double>> fragRatios, List<int> ratios,
                                        List<int> fragments)
        {
            int error = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                int sum = 0;
                for (int j = 0; j < ratios.Count; j++)
                    sum += (int)(fragRatios[j][i] * ratios[j]);
                if (sum < fragments[i])
                    error += fragments[i] - sum;
            }
            return error;
        }

        public static int FindLocalMaximumFlow(List<double> fragRatio, List<int> fragments, int sumOfIntensities)
        {
            //int bestError = int.MaxValue;
            //int bestCumul = 0;
            int cumul = 0;
            //int error = int.MaxValue;
            while (cumul <= sumOfIntensities)
            {
                //error = 0;
                for (int i = 0; i < fragments.Count; i++)
                {
                    if (fragRatio[i] * cumul > fragments[i])
                        return cumul;
                    //error += Math.Abs((int)(fragRatio[i] * cumul) - fragments[i]);
                }
                cumul++;
                /*if (error <= bestError)
                {
                    bestCumul = cumul;
                    bestError = error;
                }//*/
            }
            return cumul;// bestCumul;
        }
        
        public static int MaxFlowbkp(  List<List<double>> fragRatios, 
                                    List<int>    fragments,
                                ref List<List<double>> optimalSolutions)
        {           
            //Lists must have same number of fragments, ordered in the same manner
            int sumOfIntensities = 0;
            for (int i = 0; i < fragments.Count; i++)
                sumOfIntensities += fragments[i];                     
            
            List<int> localFlows = new List<int>();
            foreach(List<double> fragmentRatio in fragRatios)
                localFlows.Add(0);
            
            int smallestError = int.MaxValue;
            int error = Overflow(fragRatios, localFlows, fragments);
            int bestIndex = 0;
            while (bestIndex >= 0)//error <= smallestError)
            {
                bestIndex = -1;
                smallestError = error;
                for (int i = 0; i < fragRatios.Count; i++)
                {
                    localFlows[i]--;
                    int tmpError = Overflow(fragRatios, localFlows, fragments);
                    localFlows[i]++;

                    if (tmpError <= error)
                    {
                        error = tmpError;
                        bestIndex = i;
                    }
                }
                if (bestIndex >= 0)
                    localFlows[bestIndex]--;
            }
            Console.WriteLine("   Max Flow found!");
            optimalSolutions.Clear();

            int sumLocalFlows = 0;
            foreach (int localFlow in localFlows)
                sumLocalFlows += localFlow;

            List<double> newList = new List<double>();
            foreach(int localFlow in localFlows)
                newList.Add( localFlow / (double) sumLocalFlows);
            optimalSolutions.Add(newList);

            return smallestError;
        }
        
        public static int MaxFlow(List<List<double>> fragRatios,
                                    List<int> fragments,
                                ref List<List<double>> optimalSolutions)
        {
            //Lists must have same number of fragments, ordered in the same manner
            int sumOfIntensities = 0;
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

            //int smallestOverError = int.MaxValue;
            int overError = Overflow(fragRatios, localFlows, fragments);
            int underError = Underflow(fragRatios, localFlows, fragments);
            int bestIndex = 0;

            int iterSize = 1;
            while (overError > 0)//bestIndex >= 0)//error <= smallestError)//bestIndex >= 0)//
            {
                bestIndex = -1;
                int smallestUnderError = int.MaxValue;// Underflow(fragRatios, localFlows, fragments);
                for (int i = 0; i < fragRatios.Count; i++)
                //for (int i = fragRatios.Count - 1; i >= 0; i--)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= iterSize;
                        int tmpErrorMinus = Underflow(fragRatios, localFlows, fragments);

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
                /*
                if (bestIndex == -1)
                {
                    for (int i = 0; i < fragRatios.Count; i++)
                    {
                        if (localFlows[i] > 0)
                        {
                            localFlows[i]--;
                            int tmpErrorMinus = Underflow(fragRatios, localFlows, fragments);

                            if (tmpErrorMinus < error)
                            {
                                error = tmpErrorMinus;
                                bestIndex = i;
                            }
                            else
                                localFlows[i]++;
                        }
                    }
                }
                //if (bestIndex >= 0)
                //    localFlows[bestIndex]--;
                /*
                else
                {
                    for (int k = 0; k < fragRatios.Count; k++)
                    {
                        localFlows[k]--;
                        for (int i = k + 1; i < fragRatios.Count; i++)
                        {
                            localFlows[i]++;
                            int tmpErrorPlus = Overflow(fragRatios, localFlows, fragments);
                            localFlows[i]--;

                            if (tmpErrorPlus < error)
                            {
                                error = tmpErrorPlus;
                                bestIndex = i;
                            }
                        }
                        if (bestIndex >= 0)
                            localFlows[bestIndex]++;
                        else
                            localFlows[k]++;
                    }
                }//*/
                overError = Overflow(fragRatios, localFlows, fragments);
            }
            Console.WriteLine("   Max Flow found!");
            optimalSolutions.Clear();

            int sumLocalFlows = 0;
            foreach (int localFlow in localFlows)
                sumLocalFlows += localFlow;

            List<double> newList = new List<double>();
            foreach (int localFlow in localFlows)
                newList.Add(localFlow);// / (double)sumLocalFlows);
            optimalSolutions.Add(newList);

            return Underflow(fragRatios, localFlows, fragments);
        }//*/

        public static int MaxFlowBruteForce(List<double> frag2Ratio,
                                List<double> frag3Ratio,
                                List<double> frag4Ratio,
                                List<double> frag5Ratio,
                                List<int> fragments,
                                ref List<List<double>> optimalSolutions)
        {
            //Lists must have same number of fragments, ordered in the same manner
            int stepSize = 10;// 100;
            int sum = 0;
            int sumOfIntensities = 0;
            for (int i = 0; i < fragments.Count; i++)
            {
                sumOfIntensities += fragments[i];
                fragments[i] = fragments[i] / stepSize;
                sum += fragments[i];
            }
            sum += stepSize;
            int smallestError = int.MaxValue;
            List<int[]> ratios = new List<int[]>();
            int[] outputs = new int[fragments.Count];
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
                            int error = CreateFictiousSpectrum(frag2Ratio, intensity2,
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
        /*
        public static void TestMaxFlow()
        {
            Result rez = Launch();
            ProductSpectrum spectrum = rez.queries[0].spectrum;
            List<int> intensities = new List<int>();
            for (int i = 0; i < 8; i++)
                intensities.Add((int)spectrum.Peaks[i].Intensity);

            //Build fictitius ratios

            Random rd = new Random();
            List<double> ratio2 = BuildFictitiousRatios(8, rd);
            List<double> ratio3 = BuildFictitiousRatios(8, rd);
            List<double> ratio4 = BuildFictitiousRatios(8, rd);
            List<double> ratio5 = BuildFictitiousRatios(8, rd);
            List<List<double>> optimalSolutions = new List<List<double>>();
            int error = MaxFlow(ratio2, ratio3, ratio4, ratio5, intensities, ref optimalSolutions);
            vsCSVWriter writer = new vsCSVWriter(rez.dbOptions.OutputFolder + "MaxFlow_" + error + ".csv");
            writer.AddLine("Ratio 2,Ratio 3,Ratio 4,Ratio 5");
            foreach (List<double> solution in optimalSolutions)
            {
                Console.WriteLine("Max Flow computed (error of " + error + ")");
                Console.WriteLine("     Ratio 2 -> " + solution[0] + "");
                Console.WriteLine("     Ratio 3 -> " + solution[1] + "");
                Console.WriteLine("     Ratio 4 -> " + solution[2] + "");
                Console.WriteLine("     Ratio 5 -> " + solution[3] + "");
                writer.AddLine(solution[0] + "," + solution[1] + "," + solution[2] + "," + solution[3]);
            }
            writer.writeToFile();
        }//*/

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

        public static void NormalizeList(List<double> list)
        {
            double sum = 0.0;
            foreach (double item in list)
                sum += item;
            for (int i = 0; i < list.Count; i++)
                list[i] /= sum;
        }
        /*
        public static void TestGroup2MaxFlow(string fileName, List<int> fragmentIntensities)
        {
            //Y^1[6],Y^1[7],Y^1[10],C^1[5],C^1[9],C^1[12]
            //Y^1[4],B^2[5],Z^2[5],Y^1[11],A^2[12]
            //B^1[4],A^2[7],X^2[9],
            //A^1[2],B^1[3],A^1[4],B^1[4],B^2[5],B^2[6],B^2[7],Z^1[14]

            List<double> ace5 = new List<double>(new double[] { 683.4463501,1456.480347,0,1890.54895,0,0,0,0,9005.911499,0,2224.750916,0,0,2224.750916,0,4576.825195,1643.613846,1054.626343,235.686203,683.4463501 });
            NormalizeList( ace5);
            List<double> ace8 = new List<double>(new double[] { 532.7681274,0,0,416.5317688,1852.535019,2929.513916,0,1286.056152,0,0,0,0,0,2468.266693,0,4539.20166,1855.542953,0,409.5055542,532.7681274 });
            NormalizeList( ace8);
            List<double> ace12 = new List<double>(new double[] { 575.9935913,124.3237839,0,1724.810616,1824.990662,102.0265846,110.2434425,1264.954956,9028.048676,124.3237839,2752.691498,2064.453999,396.3171997,0,2009.088135,55.33430099,63.37308121,0,75.52266693,575.9935913 });
            NormalizeList( ace12);
            List<double> ace16 = new List<double>(new double[] { 2108.536354,6152.695602,2110.242016,2721.463669,106.3130417,4819.6161,891.5953293,468.4116402,125.0858727,6574.652397,206.5512619,358.8895378,3131.95974,71.50226593,574.4330292,248.7678299,679.7988434,287.7425919,565.4915924,2108.536354 });
            NormalizeList( ace16);

            List<List<double>> optimalSolutions = new List<List<double>>();
            int error = MaxFlow(ace5, ace8, ace12, ace16, fragmentIntensities, ref optimalSolutions);
            vsCSVWriter writer = new vsCSVWriter(fileName);
            writer.AddLine("Ratio 2,Ratio 3,Ratio 4,Ratio 5");
            foreach (List<double> solution in optimalSolutions)
            {
                Console.WriteLine("Max Flow computed (error of " + error + ")");
                Console.WriteLine("     ace 5 -> " + solution[0] + "");
                Console.WriteLine("     ace 8 -> " + solution[1] + "");
                Console.WriteLine("     ace 12-> " + solution[2] + "");
                Console.WriteLine("     ace 16-> " + solution[3] + "");
                writer.AddLine(solution[0] + "," + solution[1] + "," + solution[2] + "," + solution[3]);
            }
            writer.writeToFile();
        }//*/

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

            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project__Group2_All.csv";//Group 2 (all)

            //string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP11_2013\Project__Group3_All.csv";//Group 3 (all)
            
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
                List<double> tmpList = tmp.GetFragments(psmList, peptide, charge, false);                
                items.Add(tmpList);
                /*
                foreach (Precursor precursor in tmp.matchedPrecursors)
                    if (precursor.Charge == 2 && precursor.sample == sample)
                    {
                        foreach (PeptideSpectrumMatch psm in precursor.psms)
                            if (strToAim.CompareTo(psm.Peptide.Sequence) == 0)
                            {
                                peptide = psm.Peptide;
                                charge = psm.Query.precursor.Charge;
                                psmList.Add(psm);
                            }
                    }
                tmp.ExportFragmentIntensities(psmList, peptide, charge, dbOptions.OutputFolder + peptide.Sequence + ".csv");//*/
            }
            List<int> mix = new List<int>();
            for (int i = 0; i < items[0].Count; i++)
                mix.Add((int)(
                                items[0][i] * 1 + 
                                items[1][i] * 10 + 
                                items[2][i] * 20 + 
                                items[3][i] * 30));
            /*
            foreach(List<double> list in items)
                NormalizeList(list);//*/
            List<List<double>> optimalSolutions = new List<List<double>>();
            int error = MaxFlow(items, mix, ref optimalSolutions);            
            //int error = MaxFlowBruteForce(items[0], items[1], items[2], items[3], mix, ref optimalSolutions);            
            foreach (List<double> solution in optimalSolutions)
            {
                Console.WriteLine("Max Flow computed (error of " + error + ")");
                Console.WriteLine("     ace 5 -> " + solution[0] + "");
                Console.WriteLine("     ace 8 -> " + solution[1] + "");
                Console.WriteLine("     ace 12-> " + solution[2] + "");
                Console.WriteLine("     ace 16-> " + solution[3] + "");
            }
            Console.WriteLine("Number of solutions : " + optimalSolutions.Count);
            //tmp.ExportFragments(psm);
            /*
            double nbMatchingProducts = 0;
            
            foreach (Precursor precursor in tmp.matchedPrecursors)
                foreach (PeptideSpectrumMatch psm in precursor.psms)
                    if (psm.MatchingProducts > nbMatchingProducts)
                        nbMatchingProducts = psm.MatchingProducts;

            foreach (Precursor precursor in tmp.matchedPrecursors)
                foreach (PeptideSpectrumMatch psm in precursor.psms)
                    if (psm.MatchingProducts == nbMatchingProducts)
                        tmp.ExportFragments(psm);//*/
            /*
            tmp.Export(1, "All_");
            tmp.Export(0.05, "05_");
            tmp.Export(0.02, "02_");
            */
            //tmp.Export(0.05, "05_AllFragments");
            // tmp.Export(0.01, "01_");
            //tmp.Export(double.MaxValue, "All_");
            //tmp.WriteInfoToConsole();
            /*
            Optimizer op = new Optimizer(propheus);            
            op.LaunchBestPSMOptimization(tmp);//.proteins, propheus.AllQueries);
            //*/
            //Optimizer op = new Optimizer(propheus);
            //MSSearcher.Export(dbOptions.outputFolder + "5PercentOptimized_precursors.csv", Optimizer.PrecursorOptimizer(tmp.precursors, 0.05));
            //op.LaunchBestPSMOptimization(tmp);//.proteins, propheus.AllQueries);
            //op.LaunchPrecursorScoreOptimization(tmp);//.proteins, propheus.AllQueries);
            //op.Launch(tmp.proteins, propheus.AllQueries);
            /*
            propheus.Align(tmp);

            Result tmp2 = propheus.Search(1.0, false, null, propheus.CreateQueries(propheus.AllSpectras));
            tmp2.Export(0.05, "Aligned_05_");
            tmp2.Export(double.MaxValue, "Aligned_All_");
            MSSearcher.Export(dbOptions.outputFolder + "Aligned_5PercentOptimized_precursors.csv", Optimizer.PrecursorOptimizer(tmp2.precursors, 0.05));
            tmp.WriteInfoToConsole();//*/
            return tmp;
        }
    }
}
