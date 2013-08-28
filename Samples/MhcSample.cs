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
    public class MhcSample
    {
        public static DBOptions CreateOptions(string outputDir)
        {
            //@"G:\Thibault\Olivier\MnR\Databases\BD_RefGenome_WithReverse_2012-06-20.fasta";                        
            string fastaFile = @"C:\_IRIC\DATA\MHC\human_reference_2013-26-03.fasta";// Ref genome
            //@"C:\_IRIC\DATA\Test\testMHC\MHC_Sette_Peptides_20091001.fasta";//Sette Peptides
            //C:\_IRIC\DATA\MHC\4468.fasta";//MnR\Mixed_M_2013-21-01.fasta";//MHC M One sample
            

            DBOptions dbOptions = new DBOptions(fastaFile);
            dbOptions.precursorMassTolerance = new MassTolerance(8/*8withoutisotopes*/, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(0.05/*without isotopes*/, MassToleranceUnits.Da);//0.034 is a 60 000 resolution over 2000 range in mz
            dbOptions.MaximumPeptideMass = 200000;
            dbOptions.OutputFolder = outputDir;
            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            dbOptions.DigestionEnzyme = proteases["no enzyme"];
            dbOptions.NoEnzymeSearch = true;
            dbOptions.ToleratedMissedCleavages = 20;// 2;


            GraphML_List<Modification> fixMods = new GraphML_List<Modification>();
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

            dbOptions.PSMFalseDiscoveryRate = 0.01;
            dbOptions.addFragmentMods = false;
            dbOptions.addFragmentLoss = false;
            dbOptions.fragments = new Fragments();
            //dbOptions.fragments.Add(new FragmentA());
            dbOptions.fragments.Add(new FragmentB());
            //dbOptions.fragments.Add(new FragmentC());
            //dbOptions.fragments.Add(new FragmentX());
            dbOptions.fragments.Add(new FragmentY());
            //dbOptions.fragments.Add(new FragmentZ());
            return dbOptions;
        }

        public static void Launch()//Trinity.UnitTest.SettePeptideSample.Launch()
        {
            try
            {
                string projectFile = @"C:\_IRIC\DATA\MHC\project.csv"; //Patient C
                //@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUL03_2013\ProjectFile_SETTEpep_OneRAW.csv";//Sette Peptides
                //C:\_IRIC\DATA\MHC\project.csv";//MHC M One sample 
                Samples Project = new Samples(projectFile, 0);

                DBOptions dbOptions = CreateOptions( @"C:\_IRIC\DATA\Test\testMhc\");

                //ClusterOptions clusterOptions = new ClusterOptions(Project, outputDir, 5, true, 90, true);//TODO validate its in seconds for all file types
                /*
                Propheus propheus = new Propheus(dbOptions, Project);
                propheus.PrepareForSearch();
                //foreach (Sample s in propheus.AllSpectras.Keys)
                //    propheus.AllSpectras[s].tracks.Export(dbOptions.outputFolder + vsCSV.GetFileName(s.sSDF) + "_Tracks.csv");
                //Result tmp = propheus.Search(propheus.AllQueries, 1.0, false, false, null);
                
                //Save file
                Result result = new Result();
                result.dbOptions = dbOptions;
                result.queries = propheus.AllQueries;
                result.Save();
                //*/
                
                //Load file
                Result result = Result.Import(dbOptions.OutputFolder + "State.GraphML");
                Propheus propheus = new Propheus(result.dbOptions, Project);
                propheus.Load(result);
                //*/
                Result tmp = propheus.SearchLatestVersion(result.queries, true);//, 1.0, false, false, null);                
                
                tmp.WriteInfoToCsv(true);

                tmp.Export(0.05, "05_");
                tmp.Export(0.02, "02_");
                tmp.Export(1.0, "All_");

                UnitTest.Tests.MatchAllFragments(tmp);
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
                Console.WriteLine("Error in SettePeptideSample : " + ex.Message);
                Console.WriteLine(ex.StackTrace);
            }
        }
    }
}
