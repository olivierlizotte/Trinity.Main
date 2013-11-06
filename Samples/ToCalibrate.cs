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
    public class ToCalibrate
    {
        public static void Launch(bool restrain = false)
        {
            //@"G:\Thibault\Olivier\MnR\Databases\BD_RefGenome_WithReverse_2012-06-20.fasta";                        
            //Trypsin

            string outputDir = @"C:\_IRIC\DATA\Test\testPK\";
            string fastaFile = @"C:\_IRIC\DATA\Yeast\Yeast_SwissProt.fasta";//Yeast
            //@"G:\Thibault\Olivier\MQ_vs_Morpheus\Yeast_SwissProt.fasta";//Yeast
            //@"G:\Thibault\Olivier\Databases\SProHNoIso_20130430\current\sequences_2013-05-30.fa";
            //G:\Thibault\Olivier\MnR\Databases\mini_human_reference_2013-26-03.fasta";//Yeast
            string projectFile = @"C:\_IRIC\DATA\PK\NotCalibrated\project_fr1.csv";//Yeast
            //@"G:\Thibault\Olivier\MQ_vs_Morpheus\project.csv";//Yeast
            //@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JAN22_2013\_Project_FL_Single.csv";
            //G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\JUN27_2012\MR 4Rep DS\MassSense\_Test_ProjectFile_MF3.csv";
            //G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\MAR18_2013\ProjectFile_TestForProPheus.csv";

            DBOptions dbOptions = new DBOptions(fastaFile); 
            Samples Project = new Samples(projectFile, 0, dbOptions);
            dbOptions.precursorMassTolerance = new MassTolerance(80, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(80, MassToleranceUnits.ppm);//0.034 is a 60 000 resolution over 2000 range in mz
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
            if (!restrain)
            {
                //varMods.Add(ModificationDictionary.Instance["oxidation of M"]);
                //varMods.Add(ModificationDictionary.Instance["acetylation of protein N-terminus"]);
                //varMods.Add(ModificationDictionary.Instance["phosphorylation of S"]);
                //varMods.Add(ModificationDictionary.Instance["phosphorylation of T"]);
                //varMods.Add(ModificationDictionary.Instance["phosphorylation of Y"]);//*/
                varMods.Add(ModificationDictionary.Instance["Pionylation of peptide N-terminus"]);
                dbOptions.maximumVariableModificationIsoforms = 1024;// 2 * (varMods.Count + fixMods.Count);//TODO Evaluate the viability of this parameter
            }
            else
                dbOptions.maximumVariableModificationIsoforms = 2;
            dbOptions.variableModifications = varMods;

            dbOptions.addFragmentLoss = false;
            dbOptions.addFragmentMods = false;
            dbOptions.fragments = new Fragments();
            //dbOptions.fragments.Add(new FragmentA());
            dbOptions.fragments.Add(new FragmentB());
            //dbOptions.fragments.Add(new FragmentC());
            //dbOptions.fragments.Add(new FragmentX());
            dbOptions.fragments.Add(new FragmentY());
            //dbOptions.fragments.Add(new FragmentZ());

            //ClusterOptions clusterOptions = new ClusterOptions(Project, outputDir, 5, true, 90, true);//TODO validate its in seconds for all file types
            
            Propheus propheus = new Propheus(dbOptions, Project);

            dbOptions.SaveMS1Peaks = false;
            dbOptions.SaveMSMSPeaks = true;
            dbOptions.LoadSpectraIfFound = true;
            dbOptions.NbPSMToKeep = 10;
            propheus.Preload(false, false);
            propheus.PrepareQueries();

            //First pass (used to optimize parameters and score weights)
            Result tmp = propheus.SearchLatestVersion(propheus.AllQueries, false);

            //Compute number of Moded (Pionylated peptides) at given Ms tolerance
            
            long nbPsmDecoy = 0;
            long nbPsmTarget = 0;
            double avgShift = 0;
            foreach (Query query in tmp.queries)
            {
                foreach (PeptideSpectrumMatch psm in query.psms)
                {
                    if (string.IsNullOrEmpty(psm.Peptide.VariableModificationsInString))
                    {
                        nbPsmTarget++;
                        avgShift += psm.PrecursorMzError;
                    }
                    else
                        nbPsmDecoy++;
                }
            }
            avgShift /= (double)nbPsmTarget;
            dbOptions.ConSole.WriteLine("Average shift : " + avgShift);

            //tmp.WriteInfoToCsv(true);
            //tmp.Export(0.02, "FirstPass_02_");
        }
    }
}
