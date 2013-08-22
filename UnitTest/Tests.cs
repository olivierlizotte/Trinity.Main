/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Trinity;
using System.IO;
using Proteomics.Utilities;

namespace Trinity.UnitTest
{
    public class Tests
    {
        public static DBOptions CreateOptions()
        {
            DBOptions dbOptions = new DBOptions("");
            dbOptions.precursorMassTolerance = new MassTolerance(8/*8withoutisotopes*/, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(0.034/*without isotopes*/, MassToleranceUnits.Da);//0.034 is a 60 000 resolution over 2000 range in mz
            dbOptions.MaximumPeptideMass = 20000;
            dbOptions.OutputFolder = "";
            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            dbOptions.DigestionEnzyme = proteases["no enzyme"];
            dbOptions.ToleratedMissedCleavages = 20000;

            GraphML_List<Modification> fixMods = new GraphML_List<Modification>();
            dbOptions.fixedModifications = fixMods;

            GraphML_List<Modification> varMods = new GraphML_List<Modification>();
            dbOptions.variableModifications = varMods;
            dbOptions.maximumVariableModificationIsoforms = 2 * (varMods.Count + fixMods.Count);//TODO Evaluate the viability of this parameter
            dbOptions.MinimumPeptideLength = 4;
            return dbOptions;
        }

        public static bool Run()
        {
            int failed = 0;
            bool ok = true;
                        
            ok = Run_ProteinDigestTest(); Console.WriteLine("Protein Digestion test (no-enzyme) : " + (ok ? "OK" : "FAILED"));              failed += ok ? 0 : 1;
            ok = NoEnzymeDigestUnitTest.Run(); Console.WriteLine("Protein Digestion test (no-enzyme) : " + (ok ? "OK" : "FAILED")); failed += ok ? 0 : 1;

            //PushRelabelMaximumFlow ps = new PushRelabelMaximumFlow();
            //ps.test();
            //FordFulkerson ff = new FordFulkerson();
            //ff.Optimize();

            return failed == 0;
        }

        private static void AddOnce(List<string> list, string item)
        {
            foreach (string subIten in list)
                if (item.CompareTo(subIten) == 0)
                    return;
            list.Add(item);
        }

        public static void DrownCatPeptidesWithAllProteins()
        {
            //vsCSV csvPeptides = new vsCSV(@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\DEC18_2012\DMatton\Clustering_186716\Identifications.csv");
            vsCSV csvPeptides = new vsCSV(@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\DEC18_2012\DMatton\Clustering_186716\Cluster_Intensity_peptides_NormP.csv");
            vsCSVWriter writer = new vsCSVWriter(@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\DEC18_2012\DMatton\Clustering_186716\ProteinsPerPeptidesFromDatabases_AllReadingFrames.csv");

            NucleicAcid.InitHash();

            FileStream protein_fasta_database1 = new FileStream(@"G:\Thibault\Olivier\Databases\DMatton\Matton_Illumina_Anthesis_WithReverse.fasta", FileMode.Open, FileAccess.Read, FileShare.Read);
            List<Protein> proteins1 = new List<Protein>(ProteinFastaReader.ReadProteins(protein_fasta_database1, false));
            Dictionary<string, List<string>> protein1AAs = new Dictionary<string, List<string>>();
            foreach (Protein prot in proteins1)
            {
                for (int shift = 0; shift < 3; shift++)
                {
                    protein1AAs.Add(prot.Description + " | Reading Frame " + shift + " | Forward", NucleicAcid.ConvertNA3ToAAs(prot.BaseSequence, shift, false));
                    protein1AAs.Add(prot.Description + " | Reading Frame " + shift + " | Backward", NucleicAcid.ConvertNA3ToAAs(prot.BaseSequence, shift, true));
                }
            }

            FileStream protein_fasta_database2 = new FileStream(@"G:\Thibault\Olivier\Databases\DMatton\mattond_20110418_WithReverse_EditedJuly2013.fasta", FileMode.Open, FileAccess.Read, FileShare.Read);
            List<Protein> proteins2 = new List<Protein>(ProteinFastaReader.ReadProteins(protein_fasta_database2, false));
            Dictionary<string, List<string>> protein2AAs = new Dictionary<string, List<string>>();
            foreach (Protein prot in proteins2)
            {
                for (int shift = 0; shift < 3; shift++)
                {
                    protein2AAs.Add(prot.Description + " | Reading Frame " + shift + " | Forward", NucleicAcid.ConvertNA3ToAAs(prot.BaseSequence, shift, false));
                    protein2AAs.Add(prot.Description + " | Reading Frame " + shift + " | Backward", NucleicAcid.ConvertNA3ToAAs(prot.BaseSequence, shift, true));
                }
            }

            writer.AddLine(csvPeptides.LINES_LIST[0]);
            Dictionary<string, List<string>> dicOfPepProt = new Dictionary<string, List<string>>();
            for (int i = 1; i < csvPeptides.LINES_LIST.Count; i++)
            {
                string[] splits = csvPeptides.LINES_LIST[i].Split(vsCSV._Generic_Separator);
                string seq = splits[4];
                //string seq = splits[13];
                /*
                string protDesc = splits[10];
                if (protein1AAs.ContainsKey(protDesc))
                    if (!protein1AAs[protDesc].Contains(seq))
                        Console.WriteLine("Should be there 1");

                if (protein2AAs.ContainsKey(protDesc))
                    if (!protein2AAs[protDesc].Contains(seq))
                        Console.WriteLine("Should be there 1");
                //*/

                StringBuilder sb = new StringBuilder();
                foreach(string key in protein1AAs.Keys)
                    foreach(string protSeq in protein1AAs[key])
                        if(protSeq.Contains(seq))
                        {
                            sb.Append(key + ";");
                            break;
                        }
                
                foreach(string key in protein2AAs.Keys)
                    foreach(string protSeq in protein2AAs[key])
                        if(protSeq.Contains(seq))
                        {
                            sb.Append(key + ";");
                            break;
                        }

                if (sb.Length == 0)
                    Console.WriteLine("Zut");
                writer.AddLine(csvPeptides.LINES_LIST[i] + "," + sb.ToString().Trim());
            }
            writer.writeToFile();
        }

        /*
         * Too long
        public static void DrownCatPeptidesWithAllProteins()
        {            
            vsCSV csvPeptides = new vsCSV(@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\DEC18_2012\DMatton\Clustering_186716\Cluster_Intensity_peptides_NormP.csv");
            vsCSVWriter writer = new vsCSVWriter(@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\DEC18_2012\DMatton\Clustering_186716\ProteinsPerPeptidesFromDatabases.csv");
            
            NucleicAcid.InitHash();

            FileStream protein_fasta_database1 = new FileStream(@"G:\Thibault\Olivier\Databases\DMatton\Matton_Illumina_Anthesis_WithReverse.fasta", FileMode.Open, FileAccess.Read, FileShare.Read);
            List<Protein> proteins1 = new List<Protein>(ProteinFastaReader.ReadProteins(protein_fasta_database1, false));
            foreach(Protein prot in proteins1)
                prot.BaseSequence = prot.BaseSequence.Replace('T','U');

            FileStream protein_fasta_database2 = new FileStream(@"G:\Thibault\Olivier\Databases\DMatton\mattond_20110418_WithReverse_EditedJuly2013.fasta", FileMode.Open, FileAccess.Read, FileShare.Read);
            List<Protein> proteins2 = new List<Protein>(ProteinFastaReader.ReadProteins(protein_fasta_database2, false));                       
            foreach (Protein prot in proteins2)
                prot.BaseSequence = prot.BaseSequence.Replace('T', 'U'); 

            writer.AddLine(csvPeptides.LINES_LIST[0]);
            Dictionary<string, List<string>> dicOfPepProt = new Dictionary<string, List<string>>();
            for (int i = 1; i < csvPeptides.LINES_LIST.Count; i++)
            {
                string[] splits = csvPeptides.LINES_LIST[i].Split(vsCSV._Generic_Separator);
                string seq = splits[4];

                StringBuilder sb = new StringBuilder();
                foreach (string alternateSeq in NucleicAcid.ConvertAAToNAs(seq))
                {
                    foreach (Protein prot in proteins1)
                        if (prot.BaseSequence.Contains(alternateSeq))
                            sb.Append(prot.Description + " ");
                    foreach (Protein prot in proteins2)
                        if (prot.BaseSequence.Contains(alternateSeq))
                            sb.Append(prot.Description + " ");
                }
                if(sb.Length == 0)
                    Console.WriteLine("Zut");
                writer.AddLine(csvPeptides.LINES_LIST[i] + "," + sb.ToString().Trim());
            }
            writer.writeToFile();
        }//*/
        /*
         * Only one protein returned by Mascot ... so the list of protein is always of length one
        public static void DrownCatPeptidesWithAllProteins()
        {
            vsCSV csvIDs = new vsCSV(@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\DEC18_2012\DMatton\Clustering_186716\AllProteinPerPeptides.csv");
            //vsCSV csvPeptides = new vsCSV(@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\DEC18_2012\DMatton\Clustering_186716\Cluster_Intensity_peptides_NormP.csv");
            vsCSVWriter writer = new vsCSVWriter(@"G:\Thibault\-=Proteomics_Raw_Data=-\ELITE\DEC18_2012\DMatton\Clustering_186716\ProteinsPerPeptides.csv");

            Dictionary<string, List<string>> dicOfPepProt = new Dictionary<string, List<string>>();
            for (int i = 1; i < csvIDs.LINES_LIST.Count; i++)
            {
                string[] splits = csvIDs.LINES_LIST[i].Split(vsCSV._Generic_Separator);
                string seq = splits[1];
                List<string> protIDs;        
                if(!dicOfPepProt.ContainsKey(seq))
                    dicOfPepProt.Add(seq, new List<string>());
                protIDs = dicOfPepProt[seq];

                AddOnce(protIDs, splits[2]);

                for (int k = 3; k < splits.Length; k++)
                {
                    if (!string.IsNullOrEmpty(splits[k]))
                    {
                        string[] prots = splits[k].Split(' ');
                        seq = prots[3];
                        if (!dicOfPepProt.ContainsKey(seq))
                            dicOfPepProt.Add(seq, new List<string>());
                        protIDs = dicOfPepProt[seq];

                        AddOnce(protIDs, prots[0]);
                    }
                }
            }

            writer.AddLine("Peptide Sequence,Protein IDs");
            foreach (string key in dicOfPepProt.Keys)
            {
                StringBuilder sb = new StringBuilder(key);
                foreach (string prot in dicOfPepProt[key])
                    sb.Append("," + prot);
                writer.AddLine(sb.ToString());
            }
            writer.writeToFile();
        }//*/

        public static void TestGraphOptimization(Result result)
        {
            FordFulkerson ff = new FordFulkerson();
            List<FordFulkerson.Node> nodes = new List<FordFulkerson.Node>();
            List<FordFulkerson.Edge> edges = new List<FordFulkerson.Edge>();

            FordFulkerson.Node source = new FordFulkerson.Node("Source");
            nodes.Add(source);
            //foreach potential node
            //{
            //    nodes.Add(theNode);
            //    ff.AddNode(theNode);
            //    ff.AddEdge(source, theNode, capacity);
            //}
            ff.Optimize();
        }

        public static void MatchAllFragments(Result rez)
        {
            Console.WriteLine("   ************************************************************************   ");
            Console.WriteLine("   ****************** Testing All Fragments on results ********************   ");
            rez.dbOptions.fragments = new Fragments();
            rez.dbOptions.fragments.Add(new FragmentA());
            rez.dbOptions.fragments.Add(new FragmentB());
            rez.dbOptions.fragments.Add(new FragmentC());
            rez.dbOptions.fragments.Add(new FragmentX());
            rez.dbOptions.fragments.Add(new FragmentY());
            rez.dbOptions.fragments.Add(new FragmentZ());
            rez.dbOptions.addFragmentMods = true;
            rez.dbOptions.addFragmentLoss = true;
            
            foreach(Precursor precursor in rez.matchedPrecursors)
            {
                PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                PeptideSpectrumMatch psmNew = new PeptideSpectrumMatch(psm.Query, psm.Peptide, rez.dbOptions);
                precursor.psms = new PeptideSpectrumMatches();// GraphML_List<PeptideSpectrumMatch>();
                precursor.psms.Add(psmNew);
            }
            rez.WriteInfoToConsole();
            rez.Export(0.05, "05_AllFragments", true);
            Console.WriteLine("   ****************** All fragments added to best psms ********************   ");
            Console.WriteLine("   ************************************************************************   ");
        }

        public static bool Run_ProteinDigestTest()
        {
            //TODO test GPU instead
            /*
            DBOptions dbOptions = CreateOptions();
            Dictionary<string, int> sequences = new Dictionary<string,int>();
            List<Protein> proteins = Propheus.ReadProteomeFromFasta(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "UnitTest", "proteins.fasta"), false);//ETLPAMCNVYYVNCMAPLTE
            string sequence = proteins[0].BaseSequence;
            for(int i = 0; i < sequence.Length; i++)
            {
                for(int j = i + dbOptions.MinimumPeptideLength - 1; j < sequence.Length; j++)
                {
                    if(!sequences.ContainsKey(sequence.Substring(i, j - i + 1)))
                        sequences.Add(sequence.Substring(i, j - i + 1), 1);
                    else
                        sequences[sequence.Substring(i, j - i + 1)]++;
                }
            }
            foreach (Peptide peptide in Propheus.DigestProteomeOnTheFly(proteins, false, dbOptions))
                sequences[peptide.BaseSequence] -= 1;

            foreach (int val in sequences.Values)
                if (val != 0)
                    return false;//*/
            return true;
        }

        public static bool MascotCompare()
        {
            List<string> listMascotFiles = new List<string>();
            listMascotFiles.Add(@"N:\Thibault\Frederic Lamoliatte\Olivier\no MSed.csv");
            listMascotFiles.Add(@"N:\Thibault\Frederic Lamoliatte\Olivier\MSed1.csv");
            listMascotFiles.Add(@"N:\Thibault\Frederic Lamoliatte\Olivier\MSed2.csv");

            Dictionary<string, string[]> dicOfPep = new Dictionary<string, string[]>();
            for(int i = 0; i < listMascotFiles.Count; i++)
            {
                vsCSV csv = new vsCSV(listMascotFiles[i]);
                foreach (string line in csv.LINES_LIST)
                {
                    string[] splits = line.Split(vsCSV._Generic_Separator);
                    if (splits.Length > 26)
                    {
                        string key = splits[1] + "," + splits[33];
                        if (!dicOfPep.ContainsKey(key))//raw+scan+seq+mod
                            dicOfPep.Add(key, new string[3]);
                        
                        dicOfPep[key][i] = "," + splits[13] + "," + splits[14] + "," + splits[18] + "," + splits[26];
                    }
                }
            }

            vsCSVWriter writer = new vsCSVWriter(@"C:\_IRIC\DATA\Sumo\outputCompare.csv");
            foreach (string key in dicOfPep.Keys)
            {
                string str = key;
                for (int i = 0; i < listMascotFiles.Count; i++)
                    if (dicOfPep[key][i] != null)
                        str += dicOfPep[key][i];
                    else
                        str += ",,,,";
                writer.AddLine(str);
            }
            writer.writeToFile();

            return true;
        }

        public static List<string> ZincFingerProteins()
        {
            List<string> proteinID1 = new List<string>();
            List<string> proteinID2 = new List<string>();
            proteinID1.Add("P33993"); proteinID2.Add("MCM7");
            proteinID1.Add("P11388"); proteinID2.Add("TOP2A");
            proteinID1.Add("O75925"); proteinID2.Add("PIAS1");
            proteinID1.Add("Q13263"); proteinID2.Add("TRIM28");
            proteinID1.Add("Q01081"); proteinID2.Add("U2AF1");
            proteinID1.Add("P62244"); proteinID2.Add("RPS15A");
            proteinID1.Add("Q96N38"); proteinID2.Add("ZNF714");
            proteinID1.Add("Q7L590"); proteinID2.Add("MCM10");
            proteinID1.Add("P52272"); proteinID2.Add("HNRNPM");
            proteinID1.Add("Q16777"); proteinID2.Add("HIST2H2AC");
            proteinID1.Add("Q15233"); proteinID2.Add("NONO");
            proteinID1.Add("O95391"); proteinID2.Add("SLU7");
            proteinID1.Add("P55854"); proteinID2.Add("SUMO3");
            proteinID1.Add("Q9Y3A2"); proteinID2.Add("UTP11L");
            proteinID1.Add("Q9HCG1"); proteinID2.Add("ZNF160");
            proteinID1.Add("Q14592"); proteinID2.Add("ZNF460");
            proteinID1.Add("Q96IR2"); proteinID2.Add("ZNF845");
            proteinID1.Add("Q6ZR52"); proteinID2.Add("ZNF493");
            proteinID1.Add("Q03188"); proteinID2.Add("CENPC1");
            proteinID1.Add("P68431"); proteinID2.Add("HIST1H3A");
            proteinID1.Add("Q7L2R6"); proteinID2.Add("ZNF765");
            proteinID1.Add("Q96IR2"); proteinID2.Add("ZNF845");
            proteinID1.Add("Q8N4W9"); proteinID2.Add("ZNF808");
            proteinID1.Add("P62987"); proteinID2.Add("UBA52");
            proteinID1.Add("P78347"); proteinID2.Add("GTF2I");
            proteinID1.Add("P52272"); proteinID2.Add("HNRNPM");
            proteinID1.Add("Q15233"); proteinID2.Add("NONO");
            proteinID1.Add("Q08357"); proteinID2.Add("SLC20A2");
            proteinID1.Add("Q96RL1"); proteinID2.Add("UIMC1");
            proteinID1.Add("Q5QJE6"); proteinID2.Add("DNTTIP2");
            proteinID1.Add("Q08ER8"); proteinID2.Add("ZNF543");
            proteinID1.Add("Q15776"); proteinID2.Add("ZNF192");
            proteinID1.Add("Q13347"); proteinID2.Add("EIF3I");
            proteinID1.Add("P18077"); proteinID2.Add("RPL35A");
            proteinID1.Add("Q96PV6"); proteinID2.Add("LENG8");
            proteinID1.Add("O60481"); proteinID2.Add("ZIC3");
            proteinID1.Add("Q15072"); proteinID2.Add("ZNF146");
            proteinID1.Add("P58317"); proteinID2.Add("ZNF121");
            proteinID1.Add("Q99676"); proteinID2.Add("ZNF184");
            proteinID1.Add("O95391"); proteinID2.Add("SLU7");
            proteinID1.Add("P46777"); proteinID2.Add("RPL5");
            proteinID1.Add("P62269"); proteinID2.Add("RPS18");
            proteinID1.Add("Q9HCG1"); proteinID2.Add("ZNF160");
            proteinID1.Add("A2RRD8"); proteinID2.Add("ZNF320");
            proteinID1.Add("P04818"); proteinID2.Add("TYMS");
            proteinID1.Add("Q9P2P6"); proteinID2.Add("STARD9");
            proteinID1.Add("P62263"); proteinID2.Add("RPS14");
            
            proteinID1.AddRange(proteinID2);
            return proteinID1;
        }

        public static List<string> GetZincProteinsENSP()
        {
            List<string> proteinID2 = new List<string>();
            proteinID2.Add("ENSP00000473091");
            proteinID2.Add("ENSP00000375660");
            proteinID2.Add("ENSP00000287538");
            proteinID2.Add("ENSP00000359638");
            proteinID2.Add("ENSP00000249636");
            proteinID2.Add("ENSP00000297151");
            proteinID2.Add("ENSP00000428362");
            proteinID2.Add("ENSP00000315644");
            proteinID2.Add("ENSP00000411532");
            proteinID2.Add("ENSP00000419117");
            proteinID2.Add("ENSP00000393393");
            proteinID2.Add("ENSP00000344006");
            proteinID2.Add("ENSP00000307288");
            proteinID2.Add("ENSP00000346171");
            proteinID2.Add("ENSP00000359345");
            proteinID2.Add("ENSP00000325732");
            proteinID2.Add("ENSP00000325376");
            proteinID2.Add("ENSP00000330343");
            proteinID2.Add("ENSP00000468643");
            proteinID2.Add("ENSP00000326967");
            proteinID2.Add("ENSP00000318646");
            proteinID2.Add("ENSP00000458528");
            proteinID2.Add("ENSP00000457000");
            proteinID2.Add("ENSP00000458115");
            proteinID2.Add("ENSP00000385958");
            proteinID2.Add("ENSP00000428509");
            proteinID2.Add("ENSP00000311028");
            proteinID2.Add("ENSP00000385425");
            proteinID2.Add("ENSP00000211372");
            proteinID2.Add("ENSP00000416110");
            proteinID2.Add("ENSP00000412583");
            proteinID2.Add("ENSP00000403175");
            proteinID2.Add("ENSP00000393241");
            proteinID2.Add("ENSP00000472264");
            proteinID2.Add("ENSP00000396910");
            proteinID2.Add("ENSP00000471062");
            proteinID2.Add("ENSP00000388107");
            proteinID2.Add("ENSP00000470419");
            proteinID2.Add("ENSP00000471622");
            proteinID2.Add("ENSP00000472545");
            proteinID2.Add("ENSP00000470507");
            proteinID2.Add("ENSP00000473048");
            proteinID2.Add("ENSP00000471464");
            proteinID2.Add("ENSP00000358160");
            proteinID2.Add("ENSP00000329554");
            proteinID2.Add("ENSP00000352252");
            proteinID2.Add("ENSP00000350275");
            proteinID2.Add("ENSP00000244661");
            proteinID2.Add("ENSP00000439493");
            proteinID2.Add("ENSP00000366999");
            proteinID2.Add("ENSP00000367062");
            proteinID2.Add("ENSP00000353581");
            proteinID2.Add("ENSP00000444823");
            proteinID2.Add("ENSP00000439660");
            proteinID2.Add("ENSP00000322542");
            proteinID2.Add("ENSP00000322671");
            proteinID2.Add("ENSP00000322599");
            proteinID2.Add("ENSP00000387651");
            proteinID2.Add("ENSP00000418705");
            proteinID2.Add("ENSP00000420672");
            proteinID2.Add("ENSP00000369629");
            proteinID2.Add("ENSP00000291552");
            proteinID2.Add("ENSP00000381205");
            proteinID2.Add("ENSP00000273853");
            proteinID2.Add("ENSP00000340465");
            proteinID2.Add("ENSP00000429754");
            proteinID2.Add("ENSP00000429712");
            proteinID2.Add("ENSP00000322545");
            proteinID2.Add("ENSP00000253024");
            proteinID2.Add("ENSP00000342232");
            proteinID2.Add("ENSP00000362688");
            proteinID2.Add("ENSP00000353491");
            proteinID2.Add("ENSP00000392095");
            proteinID2.Add("ENSP00000400391");
            proteinID2.Add("ENSP00000276079");
            proteinID2.Add("ENSP00000362963");
            proteinID2.Add("ENSP00000362947");
            proteinID2.Add("ENSP00000441364");
            proteinID2.Add("ENSP00000332750");
            proteinID2.Add("ENSP00000402948");
            proteinID2.Add("ENSP00000332194");
            proteinID2.Add("ENSP00000228929");
            proteinID2.Add("ENSP00000411010");
            proteinID2.Add("ENSP00000376110");
            proteinID2.Add("ENSP00000347691");
            proteinID2.Add("ENSP00000424395");
            proteinID2.Add("ENSP00000379689");
            proteinID2.Add("ENSP00000470468");
            proteinID2.Add("ENSP00000367986");
            proteinID2.Add("ENSP00000418268");
            proteinID2.Add("ENSP00000420522");
            proteinID2.Add("ENSP00000352846");
            proteinID2.Add("ENSP00000470005");
            proteinID2.Add("ENSP00000388311");
            proteinID2.Add("ENSP00000291770");
            proteinID2.Add("ENSP00000318374");
            proteinID2.Add("ENSP00000459887");
            proteinID2.Add("ENSP00000366434");
            proteinID2.Add("ENSP00000427480");
            proteinID2.Add("ENSP00000421926");
            proteinID2.Add("ENSP00000366636");
            proteinID2.Add("ENSP00000211936");
            proteinID2.Add("ENSP00000470961");
            proteinID2.Add("ENSP00000406201");
            proteinID2.Add("ENSP00000409597");
            proteinID2.Add("ENSP00000290607");
            proteinID2.Add("ENSP00000362105");
            return proteinID2;
        }

        public static Dictionary<string, List<int>> ReadZNF(string filename)
        {
            vsCSV csv = new vsCSV(filename);
            Dictionary<string, List<int>> dicOfSites = new Dictionary<string,List<int>>();
            for (int i = 0; i < csv.LINES_LIST.Count; i++)
            {
                string[] splits = csv.LINES_LIST[i].Split(vsCSV._Generic_Separator);
                if (!dicOfSites.ContainsKey(splits[4] + "|" + splits[3]))
                    dicOfSites.Add(splits[4] + "|" + splits[3], new List<int>());
                dicOfSites[splits[4] + "|" + splits[3]].Add(int.Parse(splits[3]));
            }
            return dicOfSites;
        }

        public static bool LysineConservation()
        {
            List<string> zincIDs = GetZincProteinsENSP();
            Dictionary<string, List<int>> dicOfSites = ReadZNF(@"C:\_IRIC\DATA\Sumo\ZNF.csv");
            string csvMatrix = @"C:\_IRIC\DATA\Sumo\matrix_human.tsv";
            string csvToAnnotate = @"C:\_IRIC\DATA\Sumo\liste des sites SUMO.csv";
            string output = @"C:\_IRIC\DATA\Sumo\outputL.csv";
            string outputConservation = @"C:\_IRIC\DATA\Sumo\outputConservation.csv";
            string outputAll = @"C:\_IRIC\DATA\Sumo\outputConservationAll_b.csv";
            vsCSV csvM = new vsCSV(csvMatrix);
            vsCSV csvA = new vsCSV(csvToAnnotate);

            Dictionary<string, int> dicOfAnnotates = new Dictionary<string, int>();
            for (int i = 1; i < csvA.LINES_LIST.Count; i++)
            {
                string[] items = csvA.LINES_LIST[i].Split(vsCSV._Generic_Separator);
                dicOfAnnotates.Add(items[0] + items[5], i);
            }

            vsCSVWriter writer = new vsCSVWriter(output);
            writer.AddLine(csvA.getFirstLine());
                        
            Dictionary<string, List<double>> dicOfAllAA = new Dictionary<string, List<double>>();
            vsCSVWriter writerConsAll = new vsCSVWriter(outputAll);
            StringBuilder sb = new StringBuilder();
            for (char aa = 'A'; aa <= 'Z'; aa++)
            {
                dicOfAllAA.Add(aa.ToString(), new List<double>());
                sb.Append(aa + ",");
            }
            sb.Append("SpecialK,ZincEveryWhere,ZincSumo");
            dicOfAllAA.Add("SpecialK", new List<double>());
            dicOfAllAA.Add("ZincEveryWhere", new List<double>());
            dicOfAllAA.Add("ZincSumo", new List<double>());
            writerConsAll.AddLine(sb.ToString());

            writer.AddLine("AminoAcid,Some Number");
            Dictionary<string, double> dicOfAA = new Dictionary<string, double>();
            Dictionary<string, int> dicOfAANumber = new Dictionary<string, int>();

            int nbK = 0;
            List<double> echantillonageNoZero = new List<double>();
            double value = 0;
            for (int j = 0; j < csvM.LINES_LIST.Count; j++)
            {
                string[] splitsJ = csvM.LINES_LIST[j].Split('\t');
                string aa = splitsJ[3];
                if(!dicOfAA.ContainsKey(aa))
                {   
                    dicOfAANumber.Add(aa, 0);
                    dicOfAA.Add(aa, 0);
                }
                value = -1;
                if (double.TryParse(splitsJ[5], out value))
                {
                    dicOfAA[aa] += value;
                    dicOfAANumber[aa]++;

                    dicOfAllAA[aa].Add(value);
                }

                //string ensbJ = splitsJ[1];
                //int positionJ = int.Parse(splitsJ[2]);
                if ("K".CompareTo(aa) == 0)
                {
                    value = -1;
                    if (double.TryParse(splitsJ[5], out value))
                    {
                        bool found = false;
                        foreach(string id in zincIDs)
                            if(id.CompareTo(splitsJ[1]) == 0)
                                found = true;
                        if (found)
                            dicOfAllAA["ZincEveryWhere"].Add(value);

                        if (dicOfSites.ContainsKey(splitsJ[1] + "|" + splitsJ[2]))
                            dicOfAllAA["ZincSumo"].Add(value);

                        echantillonageNoZero.Add(value);

                        if (dicOfAnnotates.ContainsKey(splitsJ[1] + splitsJ[2]))
                        {
                            if (!dicOfAA.ContainsKey("SpecialK"))
                            {
                                dicOfAANumber.Add("SpecialK", 0);
                                dicOfAA.Add("SpecialK", 0);
                            }
                            dicOfAA["SpecialK"] += value;
                            dicOfAANumber["SpecialK"]++;
                            dicOfAllAA["SpecialK"].Add(value);
                            nbK++;
                            writer.AddLine(csvA.LINES_LIST[dicOfAnnotates[splitsJ[1] + splitsJ[2]]] + "," + csvM.LINES_LIST[j].Replace('\t', ','));
                        }
                    }
                }
            }
            writer.writeToFile();

            vsCSVWriter writerCons = new vsCSVWriter(outputConservation);
            foreach(string key in dicOfAA.Keys)
                writerCons.AddLine(key + "," + dicOfAA[key] / (double) dicOfAANumber[key]);
            
            double meanNoZero = 0;
            Random r = new Random();
            for (int i = 0; i < nbK; i++)
            {
                int index = (int)Math.Floor(r.NextDouble() * (echantillonageNoZero.Count - 1));
                meanNoZero += echantillonageNoZero[index];
            }
            writerCons.AddLine("Echantillonnage K," + meanNoZero / (double)nbK);
            writerCons.writeToFile();

            int lineIndex = 0;
            bool keepGoing = true;
            while (keepGoing)
            {
                StringBuilder sb2 = new StringBuilder();
                keepGoing = false;
                foreach (string key in dicOfAllAA.Keys)
                {
                    if (lineIndex < dicOfAllAA[key].Count)
                    {
                        sb2.Append(dicOfAllAA[key][lineIndex] + ",");
                        keepGoing = true;
                    }
                    else
                        sb2.Append(",");
                }
                writerConsAll.AddLine(sb2.ToString());
                lineIndex++;
            }
            writerConsAll.writeToFile();
            return true;
        }
    }
}
