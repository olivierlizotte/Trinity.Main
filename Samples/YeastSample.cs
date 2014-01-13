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
    public class YeastSample
    {
        public static void Launch(IConSol console)
        {
            //@"G:\Thibault\Olivier\MnR\Databases\BD_RefGenome_WithReverse_2012-06-20.fasta";                        
            //Trypsin

            string outputDir = @"C:\_IRIC\DATA\Yeast\Results\";
            string fastaFile = @"C:\_IRIC\DATA\Yeast\Yeast_SwissProt.fasta";//Yeast
            //@"G:\Thibault\Olivier\MQ_vs_Morpheus\Yeast_SwissProt.fasta";//Yeast
            //@"G:\Thibault\Olivier\Databases\SProHNoIso_20130430\current\sequences_2013-05-30.fa";
            //G:\Thibault\Olivier\MnR\Databases\mini_human_reference_2013-26-03.fasta";//Yeast
            string projectFile = @"C:\_IRIC\DATA\Yeast\project.csv";//Yeast
            //@"G:\Thibault\Olivier\MQ_vs_Morpheus\project.csv";//Yeast
            //@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JAN22_2013\_Project_FL_Single.csv";
            //G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUN27_2012\MR 4Rep DS\MassSense\_Test_ProjectFile_MF3.csv";
            //G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\MAR18_2013\ProjectFile_TestForProPheus.csv";
            DBOptions dbOptions = new DBOptions(fastaFile, console);
            Samples Project = new Samples(projectFile, 0, dbOptions);
            dbOptions.precursorMassTolerance = new MassTolerance(8/*8*//*8withoutisotopes*/, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(0.034/*0.034*//*without isotopes*/, MassToleranceUnits.Da);//0.034 is a 60 000 resolution over 2000 range in mz
            //dbOptions.productMassTolerance = new MassTolerance(20, MassToleranceUnits.ppm);
            dbOptions.MaximumPeptideMass = 200000;
            dbOptions.OutputFolder = outputDir;
            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            dbOptions.DigestionEnzyme = proteases["trypsin (no proline rule)"];//"no enzyme"];
            dbOptions.NoEnzymeSearch = false;// true;
            dbOptions.DecoyFusion = false;

            //dbOptions.protease = proteases["trypsin (no proline rule)"];
            dbOptions.ToleratedMissedCleavages = 2;
            dbOptions.MinimumPeptideLength = 5;
            dbOptions.MaximumPeptideLength = 300;

            GraphML_List<Modification> fixMods = new GraphML_List<Modification>();
            fixMods.Add(ModificationDictionary.Instance["carbamidomethylation of C"]);
            dbOptions.fixedModifications = fixMods;

            GraphML_List<Modification> varMods = new GraphML_List<Modification>();
            //Oxidation (M);Acetyl (Protein N-term);Phospho (STY)
            //Mods for Yeast
            varMods.Add(ModificationDictionary.Instance["oxidation of M"]);
            varMods.Add(ModificationDictionary.Instance["acetylation of protein N-terminus"]);
            varMods.Add(ModificationDictionary.Instance["phosphorylation of S"]);
            varMods.Add(ModificationDictionary.Instance["phosphorylation of T"]);
            varMods.Add(ModificationDictionary.Instance["phosphorylation of Y"]);//*/
            dbOptions.maximumVariableModificationIsoforms = 2 * (varMods.Count + fixMods.Count);//TODO Evaluate the viability of this parameter
            
            dbOptions.variableModifications = varMods;

            dbOptions.NbPSMToKeep = 16;

            dbOptions.addFragmentLoss = false;
            dbOptions.addFragmentMods = false;
            dbOptions.fragments = new Fragments();
            dbOptions.fragments.Add(new FragmentA());
            dbOptions.fragments.Add(new FragmentB());
            dbOptions.fragments.Add(new FragmentC());
            dbOptions.fragments.Add(new FragmentX());
            dbOptions.fragments.Add(new FragmentY());
            dbOptions.fragments.Add(new FragmentZ());
            
            dbOptions.dProduct = 0.0;
            dbOptions.dPrecursor = 0.1;// 0.12;
            dbOptions.dMatchingProductFraction = 0.8;// 0.45;
            dbOptions.dMatchingProduct = 0.0;// 0.5;
            dbOptions.dIntensityFraction = 0.1;// 45;// 0.0;//0.13;
            dbOptions.dIntensity = 0;
            dbOptions.dProtein = 0;
            dbOptions.dPeptideScore = 0.0;// 0.3;
            dbOptions.dFragmentScore = 0.0;// 0.5;

            //ClusterOptions clusterOptions = new ClusterOptions(Project, outputDir, 5, true, 90, true);//TODO validate its in seconds for all file types

            dbOptions.SaveMS1Peaks = true;
            dbOptions.SaveMSMSPeaks = true;
            dbOptions.LoadSpectraIfFound = true;
            Propheus propheus = new Propheus(dbOptions, Project);

            propheus.Preload(false, false);
            propheus.PrepareQueries();

            //To beat : 4653 (MaxQuant) Psm at 2%FDR
            //First pass (used to optimize parameters and score weights)
            Result tmp = propheus.SearchLatestVersion(propheus.AllQueries, true, false);//, 1.0, false, false, null);

            tmp.WriteInfoToCsv(true);
            tmp.Export(0.02, "FirstPass_02_");


            //Second search
            propheus.Preload(true);
            propheus.PrepareQueries();
            Result finalRez = propheus.SearchLatestVersion(propheus.AllQueries, false);//, 1.0, false, false, null);

            //tmp.Export(0.05, "05_");
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
        }
    }
}
