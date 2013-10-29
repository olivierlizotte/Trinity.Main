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


        public static void Optimize2()
        {
            //0.0289455865219887 : Return type (1) Precision (100000) nbProductsToKeep (13)
            //Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            //Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_MonoAce_19Oct.csv", 0);

            //0.305567296686846 : Return type (0) Precision (100000) nbProductsToKeep (5)
            //Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_DiAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            //Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_DiAce_19Oct.csv", 0);

            //0.324913505037886 : Return type (0) Precision (100000) nbProductsToKeep (5)
            Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_TriAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_TriAce_19Oct.csv", 0);

            DBOptions dbOptions = GetDBOptions(false, false);
            Result mixedResult = null;
            {
                Propheus propheus = new Propheus(dbOptions, ProjectMixed);
                propheus.Preload(false, false);
                propheus.PrepareQueries();

                mixedResult = propheus.SearchLatestVersion(propheus.AllQueries, false);
            }
            
            Result spikedResult = null;
            {
                DBOptions dbOptions2 = GetDBOptions(false, false);
                Propheus propheus2 = new Propheus(dbOptions2, ProjectRatios);

                propheus2.Preload(false, false);
                propheus2.PrepareQueries();

                spikedResult = propheus2.SearchLatestVersion(propheus2.AllQueries, false);
            }
            double smallestCumulError = double.MaxValue;
            int maxIterProducts = 25;
            for (int nbProductsToKeep = 3; nbProductsToKeep <= maxIterProducts; nbProductsToKeep++)
            {
                Console.Write("\r{0}%   ", ((100 * nbProductsToKeep) / maxIterProducts));
                for (int precision = 100000; precision <= 1000000; precision *= 10)
                {
                    for (int mflowReturnType = 0; mflowReturnType <= 2; mflowReturnType++) 
                    {
                        bool smoothedPrecursor = false;
                        List<List<double>> ratios = ComputeMaxFlow(dbOptions, ProjectRatios, ProjectMixed, mflowReturnType, nbProductsToKeep, smoothedPrecursor, precision, spikedResult, mixedResult);
                    
                        double cumulError = 0;
                        for (int i = 0; i < ratios.Count; i++)
                        {
                            double cumulSumOfRatios = 0;
                            foreach (double ratio in ratios[i])
                                cumulSumOfRatios += ratio;

                            foreach (double ratio in ratios[i])
                                cumulError += Math.Abs((ratio / cumulSumOfRatios) - (1.0 / (double)ratios[i].Count));
                        }
                        if (cumulError < smallestCumulError)
                        {
                            smallestCumulError = cumulError;
                            vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + "Cumul_BestRatios6.csv");
                            for (int i = 0; i < ratios.Count; i++)
                            {
                                string line = ProjectMixed[i].sSDF;
                                foreach (double ratio in ratios[i])
                                    line += "," + ratio;
                                writer.AddLine(line);
                            }
                            writer.AddLine(smallestCumulError + " : Return type (" + mflowReturnType + ") Precision (" + precision + ") nbProductsToKeep (" + nbProductsToKeep + ")");
                            writer.WriteToFile();
                        }
                    }
                }
            }            
            Console.WriteLine("Done!");
        }

        public static void MaxFlowThis()//string projectSingleInjections, string projectMixed)
        {
            //0.0732351971695529 : Return type (1) Precision (1000000) nbProductsToKeep (5)
            //0.0289455865219887 : Return type (1) Precision (100000) nbProductsToKeep (13)
            //Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            //Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Varied_19Oct.csv", 0);

            //0.305567296686846 : Return type (0) Precision (100000) nbProductsToKeep (5)
            //Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_DiAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            //Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_DiAce_Varied_19Oct.csv", 0);

            //0.222529420854448 : Return type (2) Precision (100000) nbProductsToKeep (19)
            Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_TriAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_TriAce_Varied_19Oct.csv", 0);

            DBOptions dbOptions = GetDBOptions(false, true);
            int mflowReturnType = 0;
            int nbProducts = 5;// 5;
            bool smoothedPrecursor = false;
            int precision = 100000;
            List<List<double>> ratios = ComputeMaxFlow(dbOptions, ProjectRatios, ProjectMixed, mflowReturnType, nbProducts, smoothedPrecursor, precision);

            vsCSVWriter writerCumul = new vsCSVWriter(dbOptions.OutputFolder + "CumulRatios6.csv");
            for (int i = 0; i < ProjectMixed.Count; i++)
            {
                string lineCumulRatio = ProjectMixed[i].sSDF;
                foreach (double n in ratios[i])
                    lineCumulRatio += "," + n;// / TotalElapsedTime;
                writerCumul.AddLine(lineCumulRatio);
            }
            //writerCumul.AddLine(strNorm);
            writerCumul.WriteToFile();
        }    
    
        public static List<List<double>> ComputeMaxFlow(
                                DBOptions dbOptions,
                                Samples ProjectRatios,
                                Samples ProjectMixed,
                                int mflowReturnType = 2,
                                int nbProducts = 8,// 5,
                                bool smoothedPrecursor = false,
                                int precision = 100000,
                                Result precomputedSpiked = null,
                                Result precomputedMixed  = null)
        {
            string baseSeq = "";
            //Types of fragments
            List<List<ProductMatch>> ratios = new List<List<ProductMatch>>();
            List<double> TrackIntensity = new List<double>();
            List<double> NormalizeFactor = new List<double>();
            CreateVirtualSpectrumFromSpikedPeptides(dbOptions, ProjectRatios, ref baseSeq, false, nbProducts,
                                                    ref ratios, ref TrackIntensity, ref NormalizeFactor, smoothedPrecursor, precomputedSpiked);
            
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
            
            List<List<double>> listOfSumOfRatio = new List<List<double>>();
            foreach(Sample sample in ProjectMixed)
            {
                //Console.WriteLine("Sample " + sample.sSDF);

                Dictionary<ProductSpectrum, bool> doneSpectrum = new Dictionary<ProductSpectrum, bool>();
                Dictionary<int, string> dicOfResults = new Dictionary<int, string>();
                List<double> sumOfRatio = new List<double>();
                foreach (double intensity in TrackIntensity)
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

                            if (keep && !DicOfSpectrum.ContainsKey(query.spectrum) && psm.Peptide.BaseSequence.CompareTo(baseSeq) == 0 && query.precursor.Charge == 2 && psm.MatchingProducts > 4)
                                DicOfSpectrum.Add(query.spectrum, psm);
                        }
                    }
                }
                List<ProductSpectrum> sortedSpectrum = new List<ProductSpectrum>(DicOfSpectrum.Keys);
                sortedSpectrum.Sort(ProductSpectrum.AscendingRetentionTimeComparison);
                double LastTimeStamp = 0;
                double TotalElapsedTime = 0;
                foreach(ProductSpectrum key in sortedSpectrum)
                {
                    PeptideSpectrumMatch psm = DicOfSpectrum[key];

                    double overFlow = 0;
                    double underFlow = 0;
                    double percentError = 0;
                    List<double> finalRatios = MaxFlowFromSpectrum(ratios, ratioNames, precision, psm.Query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                                    mflowReturnType, ref overFlow, ref underFlow, ref percentError);
                                
                    double normSumRatio = 0;
                    for (int i = 0; i < finalRatios.Count; i++)
                        normSumRatio += finalRatios[i] * NormalizeFactor[i];

                    double sumRatio = 0;
                    for (int i = 0; i < finalRatios.Count; i++)
                        sumRatio += finalRatios[i];

                    if (normSumRatio != double.NaN && normSumRatio > 0)
                    {
                        string strRatios = psm.Query.spectrum.RetentionTimeInMin.ToString() + "," + psm.Query.spectrum.PrecursorIntensity;
                        double ElapsedTime = psm.Query.spectrum.RetentionTimeInMin - LastTimeStamp;
                        for (int i = 0; i < finalRatios.Count; i++)
                        {
                            double quantifiedRatio = (finalRatios[i] / sumRatio) * psm.Query.spectrum.PrecursorIntensity * NormalizeFactor[i];// *ElapsedTime;
                            sumOfRatio[i] += quantifiedRatio;
                            strRatios += "," +quantifiedRatio;
                        }
                        for (int i = 0; i < finalRatios.Count; i++)
                            strRatios += "," + finalRatios[i];
                        for (int i = 0; i < finalRatios.Count; i++)
                            strRatios += "," + NormalizeFactor[i];

                        dicOfResults.Add(psm.Query.spectrum.ScanNumber, strRatios);
                        TotalElapsedTime += ElapsedTime;
                    }
                    LastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;
                }
                vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + sample.nameColumn + "_Ratios5.csv");
                List<int> scans = new List<int>(dicOfResults.Keys);
                scans.Sort();
                foreach (int scan in scans)
                    writer.AddLine(dicOfResults[scan]);
                writer.WriteToFile();
                listOfSumOfRatio.Add(sumOfRatio);
            }
            string strNorm = "Normalization";
            for (int i = 0; i < NormalizeFactor.Count; i++)
                strNorm += "," + NormalizeFactor[i];
            return listOfSumOfRatio;
        }

        public static void OptimizeThisMaxFlow()//string projectSingleInjections, string projectMixed)
        {
            //Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            //Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Varied_19Oct.csv", 0);

            //Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_DiAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            //Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_DiAce_Varied_19Oct.csv", 0);

            Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_TriAce_Spiked_19Oct.csv", 0);//Group 2 (all)            
            Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_TriAce_Varied_19Oct.csv", 0);

            DBOptions dbOptions = GetDBOptions(false,false);
            Propheus propheus = new Propheus(dbOptions, ProjectMixed);

            propheus.Preload(false, false);
            propheus.PrepareQueries();

            Result resultMix = propheus.SearchLatestVersion(propheus.AllQueries, false);
            
            DBOptions dbOptionsSpike = GetDBOptions(false,false);
            Propheus propheusSpike = new Propheus(dbOptionsSpike, ProjectRatios);

            propheusSpike.Preload(false, false);
            propheusSpike.PrepareQueries();

            Result precomputedSpiked = propheusSpike.SearchLatestVersion(propheusSpike.AllQueries, false);
            
            double bestSsmallestCumulError = double.MaxValue;
            List<string> bestCumulLines = new List<string>();
            for (int mflowReturnType = 0; mflowReturnType <= 2; mflowReturnType++)
                for (int precision = 100000; precision <= 100000; precision *= 10)
                    for (int considerOnlyAllThere = 1; considerOnlyAllThere <= 1; considerOnlyAllThere++)
                        for(int nbProductsToKeep = 5; nbProductsToKeep <= 15; nbProductsToKeep++)
                            for(int nbSpectras = 20; nbSpectras <= 20; nbSpectras++)
                    {
                string baseSeq = "";
                List<List<ProductMatch>> ratios = new List<List<ProductMatch>>();
                List<double> TrackIntensity = new List<double>();
                List<double> NormalizeFactor = new List<double>();
                CreateVirtualSpectrumFromSpikedPeptides(dbOptions, ProjectRatios, ref baseSeq, false, nbProductsToKeep, ref ratios, 
                                                    ref TrackIntensity, ref NormalizeFactor, false, precomputedSpiked);

                List<string> ratioNames = new List<string>();
                foreach (Sample sample in ProjectRatios)
                    ratioNames.Add(sample.nameColumn);

                //Get PSMs
                double smallestCumulError = 0.0;
                List<string> cumulLines = new List<string>();
                List<List<double>> listOfSumOfRatio = new List<List<double>>();
                foreach (Sample sample in ProjectMixed)
                {
                    string lineCumulRatio = sample.sSDF;
                    Dictionary<ProductSpectrum, bool> doneSpectrum = new Dictionary<ProductSpectrum, bool>();
                    Dictionary<int, string> dicOfResults = new Dictionary<int, string>();
                    List<double> sumOfRatio = new List<double>();
                    foreach (double intensity in TrackIntensity)
                        sumOfRatio.Add(0);

                    foreach (Query query in resultMix.queries)
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

                                    double overFlow = 0;
                                    double underFlow = 0;
                                    double percentError = 0;
                                    List<double> finalRatios = MaxFlowFromSpectrum(ratios, ratioNames, precision, query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                               mflowReturnType, ref overFlow, ref underFlow, ref percentError);

                                    double normSumRatio = 0;
                                    for (int i = 0; i < finalRatios.Count; i++)
                                        normSumRatio += finalRatios[i] * NormalizeFactor[i];
                                    double sumRatio = 0;
                                    for (int i = 0; i < finalRatios.Count; i++)
                                        sumRatio += finalRatios[i];

                                    if (!(double.IsNaN(normSumRatio) || normSumRatio == 0))
                                    {
                                        if (considerOnlyAllThere == 0)
                                        {
                                            bool AllThere = true;
                                            foreach (double dbl in finalRatios)
                                                if (dbl == 0)
                                                    AllThere = false;
                                            if (AllThere)
                                            {
                                                for (int i = 0; i < finalRatios.Count; i++)
                                                    sumOfRatio[i] += (finalRatios[i] / sumRatio) * query.spectrum.PrecursorIntensity * NormalizeFactor[i];
                                                    //sumOfRatio[i] += ((finalRatios[i] * NormalizeFactor[i]) / normSumRatio) * query.spectrum.PrecursorIntensity;
                                            }
                                        }
                                        else
                                        {
                                            for (int i = 0; i < finalRatios.Count; i++)
                                                sumOfRatio[i] += (finalRatios[i] / sumRatio) * query.spectrum.PrecursorIntensity * NormalizeFactor[i];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /*
0.058823529	0.941176471
0.111111111	0.888888889
0.2	0.8
0.333333333	0.666666667
0.5	0.5
0.666666667	0.333333333
0.8	0.2
*/

                    double cumulSumOfRatios = 0;
                    foreach (double ratio in sumOfRatio)
                        cumulSumOfRatios += ratio;
                    if (cumulSumOfRatios == 0)
                        cumulSumOfRatios = 0.5;
                    else
                        cumulSumOfRatios *= 0.5;
                    
                    for (int indexPeptide = 0; indexPeptide < sumOfRatio.Count; indexPeptide++)
                    {
                        switch (sample.nameColumn)
                        {
                            case "5":
                                if (indexPeptide == 0 || indexPeptide == 2 || indexPeptide == 4)
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.058823529);
                                else
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.941176471);
                                break;
                            case "10":
                                if (indexPeptide == 0 || indexPeptide == 2 || indexPeptide == 4)
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.111111111);
                                else
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.888888889);
                                break;
                            case "20":
                                if (indexPeptide == 0 || indexPeptide == 2 || indexPeptide == 4)
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.2);
                                else
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.8);
                                break;
                            case "40":
                                if (indexPeptide == 0 || indexPeptide == 2 || indexPeptide == 4)
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.333333333);
                                else
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.666666667);
                                break;
                            case "80":
                                if (indexPeptide == 0 || indexPeptide == 2 || indexPeptide == 4)
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.5);
                                else
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.5);
                                break;
                            case "160":
                                if (indexPeptide == 0 || indexPeptide == 2 || indexPeptide == 4)
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.666666667);
                                else
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.333333333);
                                break;
                            case "320":
                                if (indexPeptide == 0 || indexPeptide == 2 || indexPeptide == 4)
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.8);
                                else
                                    smallestCumulError += Math.Abs((sumOfRatio[indexPeptide] / cumulSumOfRatios) - 0.2);
                                break;
                        }
                    }
                    foreach (double n in sumOfRatio)
                        lineCumulRatio += "," + n;
                    cumulLines.Add(lineCumulRatio);
                }
                if (smallestCumulError < bestSsmallestCumulError)
                {
                    bestSsmallestCumulError = smallestCumulError;
                    bestCumulLines = cumulLines;
                    bestCumulLines.Add(bestSsmallestCumulError + " : Return type (" + mflowReturnType + ") Precision (" + precision + " ConsiderOnlyAllThere (" + considerOnlyAllThere + ") nbSpectras (" + nbSpectras + ") nbProductsToKeep (" + nbProductsToKeep + ")");

                    vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + "Cumul_BestRatios.csv");
                    foreach (string line in bestCumulLines)
                        writer.AddLine(line);
                    string strNorm = "Normalization";
                    for (int i = 0; i < NormalizeFactor.Count; i++)
                        strNorm += "," + NormalizeFactor[i];
                    writer.AddLine(strNorm);
                    writer.WriteToFile();
                    Console.WriteLine(bestSsmallestCumulError + " : Return type (" + mflowReturnType + ") Precision (" + precision + " ConsiderOnlyAllThere (" + considerOnlyAllThere + ") nbSpectras (" + nbSpectras + ") nbProductsToKeep (" + nbProductsToKeep + ")");
                }
            }
            Console.WriteLine("Done!");
            /*
            vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + "Cumul_BestRatios.csv");
            foreach (string line in bestCumulLines)
                writer.AddLine(line);
            writer.writeToFile();//*/
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

        public static List<double> MaxFlowFromSpectrum(List<List<ProductMatch>> ratiosToFit, List<string> ratioNames, 
                                            int precision, GraphML_List<MsMsPeak> capacity, MassTolerance tolerance, 
                                            int returnType,//0 for max flow, 1 for best flow, 2 for average
                                            ref double overFlow, ref double underFlow, ref double errorInPercent)
        {
            List<List<double>> solutions = new List<List<double>>();
            List<long> average = new List<long>();
            List<MsMsPeak> expandedCapacity = new List<MsMsPeak>();
            foreach (MsMsPeak peak in capacity)
                expandedCapacity.Add(new MsMsPeak(peak.MZ, peak.Intensity * precision, peak.Charge));

            double error = ComputeMaxFlow(ratiosToFit, expandedCapacity, tolerance, ref solutions, ref errorInPercent, ref average);

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
                    ref List<List<ProductMatch>> FinalSpikedProducts, ref List<double> TrackIntensity, ref List<double> Normalizor, bool smoothPrecursor, Result precomputed = null)
        {   
            //Annotate spectrums
            Result tmp = precomputed;
            if (precomputed == null)
            {
                //DBOptions dbOptions = GetDBOptions(loadFromRaw, false);
                Propheus propheus = new Propheus(dbOptions, Project);

                propheus.Preload(false, false);
                propheus.PrepareQueries();

                tmp = propheus.SearchLatestVersion(propheus.AllQueries, false);
                tmp.WriteInfoToCsv(false);
            }

            int charge = 2;
            Peptide peptide = null;
            FinalSpikedProducts.Clear();
            List<List<ProductMatch>> SpikedProducts = new List<List<ProductMatch>>();

            //Get list of PSMs per sample
            double averagePrecursorArea = 0;
            List<List<PeptideSpectrumMatch>> listOfPSMs = new List<List<PeptideSpectrumMatch>>();
            foreach (Sample sample in Project)
            {
                List<PeptideSpectrumMatch> psmList = new List<PeptideSpectrumMatch>();
                foreach (Query query in tmp.queries)
                {
                    if (query.precursor.Charge == charge && query.precursor.sample == sample)
                    {
                        foreach (PeptideSpectrumMatch psm in query.psms)
                        {
                            if (sample.nameColumn.CompareTo(psm.Peptide.Sequence) == 0 && psm.MatchingProducts > 4)//.ProbabilityScore() > 0.020)// && psm.ProbabilityScore() > bestPsmScore)
                            {
                                peptide = psm.Peptide;
                                BaseSequence = peptide.BaseSequence;
                                psmList.Add(psm);
                            }
                        }
                    }
                }

                psmList.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
                listOfPSMs.Add(psmList);
                List<ProductMatch> productList = GetCombinedSpectrum(tmp.dbOptions, psmList, peptide, charge);

                productList.Sort(ProductMatch.AscendingWeightComparison);
                //if (productList.Count > nbProductsToKeep)
                //    productList.RemoveRange(0, productList.Count - nbProductsToKeep);//*/

                double PrecursorArea = ComputePrecursorArea(psmList, smoothPrecursor);
                averagePrecursorArea += PrecursorArea;

                SpikedProducts.Add(productList);
                TrackIntensity.Add(PrecursorArea);// / (double)nbPrec);
                //tmp.ExportFragmentIntensitiesForAllPSM(psmList, peptide, charge, tmp.dbOptions.OutputFolder + sample.nameColumn + ".csv");
            }
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
                FinalSpikedProducts.Add(GetCombinedSpectrum(tmp.dbOptions, listOfPSMs[i], peptide, charge, DicOfFragmentsToKeep));            
            
            for (int i = 0; i < Project.Count; i++)
            {
                double factor = averagePrecursorArea / TrackIntensity[i];// ((averagePrecursorIntensity / TrackIntensity[i]) - 1.0) * 0.1 + 1.0;
                for (int j = 0; j < FinalSpikedProducts[i].Count; j++)
                    FinalSpikedProducts[i][j].obsIntensity *= factor;//*/
                Normalizor.Add(factor);
            }
        }

        public static double ComputePrecursorArea(List<PeptideSpectrumMatch> psms, bool smooth)
        {
            double lastTimeStamp = 0;
            List<double> timeGap = new List<double>();
            List<double> precursorIntensities = new List<double>();
            foreach (PeptideSpectrumMatch psm in psms)
            {
                if (psm.Query.spectrum.PrecursorIntensity > 0 && lastTimeStamp > 0)
                {
                    timeGap.Add(psm.Query.spectrum.RetentionTimeInMin - lastTimeStamp);
                    precursorIntensities.Add(psm.Query.spectrum.PrecursorIntensity);
                }
                lastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;
            }
            //Smooth the curve

            if (smooth)
            {
                List<double> timeGapSmooth = new List<double>();
                List<double> precursorIntensitiesSmoothed = new List<double>();
                for (int i = 0; i < timeGap.Count; i++)
                {
                    double cumulIntensity = 0;
                    int nbItems = 0;
                    for (int k = (i - 2 < 0 ? 0 : i - 2); k < timeGap.Count - 2; k++)
                    {
                        nbItems++;
                        cumulIntensity += precursorIntensities[k];
                    }
                    if (nbItems > 0)
                    {
                        timeGapSmooth.Add(timeGap[0]);
                        precursorIntensitiesSmoothed.Add(cumulIntensity / (double)nbItems);
                    }
                }

                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGapSmooth.Count; i++)
                    fragSpectrumArea += timeGapSmooth[i] * precursorIntensitiesSmoothed[i];
                return fragSpectrumArea;
            }
            else
            {
                double fragSpectrumArea = 0;
                for (int i = 0; i < timeGap.Count; i++)
                    fragSpectrumArea += timeGap[i] * precursorIntensities[i];
                return fragSpectrumArea;
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
                foreach (ProductMatch pm in matchList)
                    newPM.obsIntensity += pm.obsIntensity + pm.obsIntensity * DicOfProductMsMsFactor[pm];

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
                        newMatch.weight *= newMatch.obsIntensity;
                        products.Add(newMatch);
                    }
                }
            }
            return products;
        }
    }
}
