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
    public class SettePeptideSample
    {
        public static void Launch()//Trinity.UnitTest.SettePeptideSample.Launch()
        {
            string outputDir = @"C:\_IRIC\DATA\Test\testMHCSette\";
            string fastaFile = @"C:\_IRIC\DATA\MHC Sette\MHC_Sette_Peptides_20091001.fasta";
            string projectFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUL03_2013\ProjectFile_SETTEpep_OneRAW.csv";

            DBOptions dbOptions = new DBOptions(fastaFile);
            try
            {       
                Samples Project = new Samples(projectFile, 0, dbOptions);
                dbOptions.precursorMassTolerance = new MassTolerance(5, MassToleranceUnits.ppm);
                dbOptions.productMassTolerance = new MassTolerance(0.068, MassToleranceUnits.Da);//0.034 is a 60 000 resolution over 2000 range in mz
                dbOptions.MaximumPeptideMass = 200000;
                dbOptions.OutputFolder = outputDir;
                ProteaseDictionary proteases = ProteaseDictionary.Instance;
                dbOptions.DigestionEnzyme = proteases["no enzyme"];
                dbOptions.NoEnzymeSearch = false;
                dbOptions.ToleratedMissedCleavages = 20;// 2;


                GraphML_List<Modification> fixMods = new GraphML_List<Modification>();
                    //fixMods.Add(ModificationDictionary.Instance["carbamidomethylation of C"]);
                dbOptions.fixedModifications = fixMods;

                GraphML_List<Modification> varMods = new GraphML_List<Modification>();
                    varMods.Add(ModificationDictionary.Instance["oxidation of M"]);//+ Deamidation M Q
                    varMods.Add(ModificationDictionary.Instance["phosphorylation of S"]);
                    varMods.Add(ModificationDictionary.Instance["phosphorylation of T"]);
                    varMods.Add(ModificationDictionary.Instance["phosphorylation of Y"]);//*/
                    varMods.Add(ModificationDictionary.Instance["deamidation of N"]);
                    varMods.Add(ModificationDictionary.Instance["deamidation of Q"]);
                    varMods.Add(ModificationDictionary.Instance["cysteinylation of C"]);
                dbOptions.variableModifications = varMods;

                dbOptions.maximumVariableModificationIsoforms = 1024;// 2 * (varMods.Count + fixMods.Count);//TODO Evaluate the viability of this parameter

                dbOptions.PSMFalseDiscoveryRate = 0.05;
                dbOptions.addFragmentMods = false;
                dbOptions.addFragmentLoss = false;
                dbOptions.fragments = new Fragments();
//                dbOptions.fragments.Add(new FragmentA());
                dbOptions.fragments.Add(new FragmentB());
//                dbOptions.fragments.Add(new FragmentC());
//                dbOptions.fragments.Add(new FragmentX());
                dbOptions.fragments.Add(new FragmentY());
//                dbOptions.fragments.Add(new FragmentZ());

                //ClusterOptions clusterOptions = new ClusterOptions(Project, outputDir, 5, true, 90, true);//TODO validate its in seconds for all file types

                Propheus propheus = new Propheus(dbOptions, Project);

                dbOptions.SaveMS1Peaks = false;
                dbOptions.SaveMSMSPeaks = false;
                dbOptions.LoadSpectraIfFound = true;
                propheus.Preload(true);
                propheus.PrepareQueries();
                
                //First pass (used to optimize parameters and score weights)
                Result tmp = propheus.SearchVersionAugust2013(propheus.AllQueries, true);//, 1.0, false, false, null);

                tmp.WriteInfoToCsv(true);
                tmp.Export(0.02, "FirstPass_02_");

                dbOptions.SaveMS1Peaks = true;
                dbOptions.SaveMSMSPeaks = true;
                
                //Second search
                propheus.Preload(true);
                propheus.PrepareQueries();
                Result finalRez = propheus.SearchLatestVersion(propheus.AllQueries, false);//, 1.0, false, false, null);
                
                //tmp.Export(0.05, "05_");
                tmp.Export(0.02, "02_");
                //tmp.Export(1.0, "All_");

                //UnitTest.Tests.MatchAllFragments(tmp);
                //tmp.WriteInfoToConsole();
                //tmp.Save();
                /*
                Optimizer op = new Optimizer(propheus);
                op.LaunchBestPSMOptimization(tmp);
                op.LaunchPrecursorScoreOptimization(tmp);
                /*
                propheus.Align(tmp);

                Result tmp2 = propheus.Search(1.0, false, null, propheus.CreateQueries(propheus.AllSpectras));
                tmp2.Export(0.05, "Aligned_05_");
                tmp2.Export(0.02, "Aligned_02_");
                tmp2.Export(double.MaxValue, "Aligned_All_");
                MSSearcher.Export(dbOptions.outputFolder + "Aligned_5PercentOptimized_precursors.csv", Optimizer.PrecursorOptimizer(tmp2.precursors, 0.05));
                tmp2.WriteInfoToConsole();//*/
            }
            catch (Exception ex)
            {
                dbOptions.ConSole.WriteLine("Error in SettePeptideSample : " + ex.Message);
                dbOptions.ConSole.WriteLine(ex.StackTrace);
            }
        }
    }
}
