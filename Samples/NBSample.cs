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
        public static void CreateFictiousSpectrum(  List<double> frag2Ratio, int out2Ratio,
                                        List<double> frag3Ratio, int out3Ratio,
                                        List<double> frag4Ratio, int out4Ratio,
                                        List<double> frag5Ratio, int out5Ratio,
                                        ref int[] output)
        {
            for (int i = 0; i < frag2Ratio.Count; i++)
            {
                int sum = 0;
                sum += (int)(frag2Ratio[i] * out2Ratio);
                sum += (int)(frag3Ratio[i] * out3Ratio);
                sum += (int)(frag4Ratio[i] * out4Ratio);
                sum += (int)(frag5Ratio[i] * out5Ratio);
                output[i] = sum;
            }
        }

        public static int FindMaximumFlow(List<double> fragRatio, List<int> fragments, int sumOfIntensities)
        {
            int bestError = int.MaxValue;            
            int cumul = 0;
            int error = int.MaxValue;
            while (error <= bestError && cumul < sumOfIntensities)
            {
                bestError = error;
                error = 0;
                for (int i = 0; i < fragments.Count; i++)
                    error += Math.Abs((int)(fragRatio[i] * cumul) - fragments[i]);
                cumul++;
            }
            return cumul;
        }

        public static int MaxFlow(List<double> frag2Ratio, 
                                List<double> frag3Ratio,
                                List<double> frag4Ratio,
                                List<double> frag5Ratio,
                                List<int>    fragments,
                                ref List<List<double>> optimalSolutions)
        {
            //Lists must have same number of fragments, ordered in the same manner
            int stepSize = 20;
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
            List<List<int>> ratios = new List<List<int>>();
            int[] outputs = new int[fragments.Count];
            int localMax2 = FindMaximumFlow(frag2Ratio, fragments, sum);
            int localMax3 = FindMaximumFlow(frag3Ratio, fragments, sum);
            int localMax4 = FindMaximumFlow(frag4Ratio, fragments, sum);
            int localMax5 = FindMaximumFlow(frag5Ratio, fragments, sum);
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
                            CreateFictiousSpectrum( frag2Ratio, intensity2,
                                                    frag3Ratio, intensity3,
                                                    frag4Ratio, intensity4,
                                                    frag5Ratio, intensity5,
                                                    ref outputs);
                            int error = 0;
                            for (int i = 0; i < outputs.Length; i++)
                                error += Math.Abs(outputs[i] - fragments[i]);
                            if (error <= smallestError)
                            {
                                bool isOk = true;
                                if (error == 0)
                                    isOk = true;
                                if (isOk)
                                {
                                    List<int> newRatios = new List<int>();
                                    newRatios.Add(intensity2);
                                    newRatios.Add(intensity3);
                                    newRatios.Add(intensity4);
                                    newRatios.Add(intensity5);
                                    if (error == smallestError)
                                        ratios.Add(newRatios);
                                    else
                                    {
                                        ratios.Clear();
                                        ratios.Add(newRatios);
                                    }
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
            for(int i = 0; i < ratios.Count; i++)
            {
                List<double> newList = new List<double>();
                for(int j = 0; j < ratios[i].Count; j++)
                    newList.Add( ratios[i][j] / (double) sumOfIntensities);
                optimalSolutions.Add(newList);
            }

            return smallestError * stepSize;// + sum;
        }

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

        public static Result Launch(bool restrain = false)
        {
            //@"G:\Thibault\Olivier\MnR\Databases\BD_RefGenome_WithReverse_2012-06-20.fasta";                        
            //Trypsin

            string outputDir = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\FEB13_2013\MSMS files\Trinity_Output\";
            string fastaFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\FEB13_2013\MSMS files\Peptide.fasta";//Yeast
            //@"G:\Thibault\Olivier\MQ_vs_Morpheus\Yeast_SwissProt.fasta";//Yeast
            //@"G:\Thibault\Olivier\Databases\SProHNoIso_20130430\current\sequences_2013-05-30.fa";
            //G:\Thibault\Olivier\MnR\Databases\mini_human_reference_2013-26-03.fasta";//Yeast
            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\FEB13_2013\MSMS files\Project_2.csv";//Yeast
            //@"G:\Thibault\Olivier\MQ_vs_Morpheus\project.csv";//Yeast
            //@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JAN22_2013\_Project_FL_Single.csv";
            //G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUN27_2012\MR 4Rep DS\MassSense\_Test_ProjectFile_MF3.csv";
            //G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\MAR18_2013\ProjectFile_TestForProPheus.csv";
            Samples Project = new Samples(projectFile, 0);
            DBOptions dbOptions = new DBOptions(fastaFile);
            dbOptions.precursorMassTolerance = new MassTolerance(30/*8*//*8withoutisotopes*/, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(0.05/*0.034*//*without isotopes*/, MassToleranceUnits.Da);//0.034 is a 60 000 resolution over 2000 range in mz
            dbOptions.MaximumPeptideMass = 200000;
            dbOptions.OutputFolder = outputDir;
            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            dbOptions.DigestionEnzyme = proteases["no enzyme"];//trypsin (no proline rule)"];
            dbOptions.NoEnzymeSearch = false;// true;
            dbOptions.DecoyFusion = false;

            //dbOptions.protease = proteases["trypsin (no proline rule)"];
            dbOptions.ToleratedMissedCleavages = 200;// 2;
            dbOptions.MinimumPeptideLength = 5;
            dbOptions.MaximumPeptideLength = 300;

            GraphML_List<Modification> fixMods = new GraphML_List<Modification>();
            //fixMods.Add(ModificationDictionary.Instance["propionylation of K"]);
            dbOptions.fixedModifications = fixMods;

            GraphML_List<Modification> varMods = new GraphML_List<Modification>();
            //Oxidation (M);Acetyl (Protein N-term);Phospho (STY)
            //Mods for Yeast
            if (!restrain)
            {
                varMods.Add(ModificationDictionary.Instance["acetylation of K"]);
                varMods.Add(ModificationDictionary.Instance["propionylation of K"]);
                dbOptions.maximumVariableModificationIsoforms = 1024;// 2 * (varMods.Count + fixMods.Count);//TODO Evaluate the viability of this parameter
            }
            else
                dbOptions.maximumVariableModificationIsoforms = 2;
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
            
            Propheus propheus = new Propheus(dbOptions, Project);

            propheus.PrepareForSearch();

            Result tmp = propheus.SearchLatestVersion(propheus.AllQueries, true);
            //Result tmp = propheus.Search(propheus.AllQueries, 1.0, false, false, null);
            //tmp.Export(0.05, "05_");

            //tmp.Save();//*/
            //Result tmp = Result.Import(dbOptions.outputFolder + "State.GraphML");
            

            //MSSearcher.Export(dbOptions.outputFolder + "TESTOptimizedV2_precursors.csv", OptimizerV2.PrecursorOptimizer(tmp.precursors, 0.05));

            //UnitTest.Tests.MatchAllFragments(tmp); 
            tmp.WriteInfoToCsv(false);
                        
            foreach (Precursor precursor in tmp.matchedPrecursors)
                foreach (PeptideSpectrumMatch psm in precursor.psms)
                    if ("GKGGKGLGKGGAKR".CompareTo(psm.Peptide.BaseSequence) == 0)
                        tmp.ExportFragments(psm);
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

            tmp.Export(1, "All_");
            tmp.Export(0.05, "05_");
            tmp.Export(0.02, "02_");
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
