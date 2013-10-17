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
    public class RT_MHC
    {
        public static void Launch()//Trinity.UnitTest.SettePeptideSample.Launch()
        {
            try
            {       
                string outputDir = @"C:\_IRIC\DATA\Test\testRTMHC\";
                string fastaFile = @"C:\_IRIC\DATA\MHC Sette\MHC_Sette_Peptides_20091001.fasta";
                string projectFile = //@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\AUG06_2013\RT_MHC\Project_TEST_600mM.csv"; 
                    @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\MAR26_2013\Project_NonFAIMS.csv";
                    //                 @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\SEP10_2013\Project_TEST_75_100_300mM.csv";
                Samples Project = new Samples(projectFile, 0);
                DBOptions dbOptions = new DBOptions(fastaFile);
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
                dbOptions.fragments.Add(new FragmentB());
                dbOptions.fragments.Add(new FragmentY());

                dbOptions.MinimumPrecursorIntensityRatioInIsolationWindow = 0.05;
                Propheus propheus = new Propheus(dbOptions, Project);
                
                dbOptions.SaveMS1Peaks = true;
                dbOptions.SaveMSMSPeaks = true;
                dbOptions.LoadSpectraIfFound = true;
                propheus.Preload(true);
                propheus.PrepareQueries();                
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error in SettePeptideSample : " + ex.Message);
                Console.WriteLine(ex.StackTrace);
            }
        }
    }
}
