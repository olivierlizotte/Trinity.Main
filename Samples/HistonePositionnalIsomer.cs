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
    public class HistonePositionnalIsomer
    {
        public static DBOptions GetDBOptions(bool loadFromRaw, bool onlyYions, IConSol console)
        {
            string outputDir = @"C:\_IRIC\DATA\Test\testNB\Iso3\";
            string fastaFile = @"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\FEB13_2013\MSMS files\Peptide.fasta";
            DBOptions dbOptions = new DBOptions(fastaFile, console);
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

        public static void MonoAce(IConSol console)
        {
            DBOptions dbOptions = GetDBOptions(false, false, console);
            Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Spiked_19Oct.csv", 0, dbOptions);//Group 2 (all)            
            Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_MonoAce_Varied_19Oct.csv", 0, dbOptions);
            Samples ProjectStable = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_MonoAce_19Oct.csv", 0, dbOptions);            
            
            int nbProductsMin = 4;// 4;
            int nbProductsMax = 8;// 14;
            bool smoothedPrecursor = false;
            int precision = 10000;
            int maxCharge = 2;

            PositionnalIsomerSolver.Solve(ProjectRatios, ProjectStable, ProjectMixed, dbOptions, nbProductsMin, nbProductsMax, smoothedPrecursor, precision, maxCharge);
        }

        public static void DiAce(IConSol console)
        {
            DBOptions dbOptions = GetDBOptions(false, false, console);
            Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_DiAce_Spiked_19Oct.csv", 0, dbOptions);
            Samples ProjectMixed =  new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_DiAce_Varied_19Oct.csv", 0, dbOptions);
            Samples ProjectStable = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_DiAce_19Oct.csv", 0, dbOptions);

            int nbProductsMin = 4;// 4;
            int nbProductsMax = 8;// 14;
            bool smoothedPrecursor = false;
            int precision = 10000;
            int maxCharge = 2;

            PositionnalIsomerSolver.Solve(ProjectRatios, ProjectStable, ProjectMixed, dbOptions, nbProductsMin, nbProductsMax, smoothedPrecursor, precision, maxCharge);
        }

        public static void TriAce(IConSol console)
        {
            DBOptions dbOptions = GetDBOptions(false, false, console);
            Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_TriAce_Spiked_19Oct.csv", 0, dbOptions);
            Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_TriAce_Varied_19Oct.csv", 0, dbOptions);
            Samples ProjectStable = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_TriAce_19Oct.csv", 0, dbOptions);

            int nbProductsMin = 4;// 4;
            int nbProductsMax = 8;// 14;
            bool smoothedPrecursor = false;
            int precision = 10000;
            int maxCharge = 2;

            PositionnalIsomerSolver.Solve(ProjectRatios, ProjectStable, ProjectMixed, dbOptions, nbProductsMin, nbProductsMax, smoothedPrecursor, precision, maxCharge);
        }

        public static void Mixed(IConSol console)
        {
            DBOptions dbOptions = GetDBOptions(false, false, console);
            Samples ProjectRatios = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_AllAce_Spiked_19Oct.csv", 0, dbOptions);
            Samples ProjectMixed = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_AllAce_Varied_19Oct.csv", 0, dbOptions);
            Samples ProjectStable = new Samples(@"C:\_IRIC\DATA\NB\ProjectTest_StableMix_AllAce_19Oct.csv", 0, dbOptions);

            int nbProductsMin = 4;// 4;
            int nbProductsMax = 8;// 14;
            bool smoothedPrecursor = false;
            int precision = 10000;
            int maxCharge = 2;

            PositionnalIsomerSolver.Solve(ProjectRatios, ProjectStable, ProjectMixed, dbOptions, nbProductsMin, nbProductsMax, smoothedPrecursor, precision, maxCharge);
        }
    }
}
