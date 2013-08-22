/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Proteomics.Utilities;

namespace Trinity.UnitTest
{
    public static class NoEnzymeDigestUnitTest
    {
        public static bool Run()
        {
            //TODO test GPU instead            
            DBOptions dbOptions = MhcSample.CreateOptions("");
            Dictionary<string, int> sequences = new Dictionary<string, int>();
            List<Protein> proteins = Propheus.ReadProteomeFromFasta(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "UnitTest", "proteins.fasta"), false);//ETLPAMCNVYYVNCMAPLTE
            string sequence = proteins[0].BaseSequence;
            double[] proteinMasses = new double[sequence.Length];

            List<double> precursors = new List<double>();
            for (int i = 0; i < sequence.Length; i++)
            {
                for (int j = i + dbOptions.MinimumPeptideLength - 1; j < sequence.Length; j++)
                {
                    int size = j - i + 1;
                    if (size <= dbOptions.MaximumPeptideLength)
                    {
                        string subStr = sequence.Substring(i, j - i + 1);
                        if (!sequences.ContainsKey(subStr))
                            sequences.Add(subStr, 1);
                        else
                            sequences[subStr]++;

                        double mass = Constants.WATER_MONOISOTOPIC_MASS;
                        for (int k = 0; k < subStr.Length; k++)
                            mass += AminoAcidMasses.GetMonoisotopicMass(subStr[k]);
                        precursors.Add(mass);
                    }
                }
                proteinMasses[i] = AminoAcidMasses.GetMonoisotopicMass(sequence[i]);
            }
            precursors.Sort();

            Queries queries = new Queries(dbOptions, precursors.ToArray());
            Digestion ps = new Digestion(dbOptions);
            List<Protein> lProt = new List<Protein>();
            lProt.Add(proteins[0]); 
            //for each protein, build matrix of mass
            //Trinity_Gpu.ProteinDigest pg = new Trinity_Gpu.ProteinDigest(precursors.ToArray(), sequence.Length);
            //Test twice to test that precursor list stays in gpu memory
            for (int iter = 0; iter < 2; iter++)
            {
                Dictionary<string, int> sequencesTmp = new Dictionary<string, int>(sequences);

                foreach (Tuple<Peptide, int> item in ps.DigestProteomeOnTheFlyNoEnzyme(lProt, queries))
                {
                    sequencesTmp[item.Item1.BaseSequence] -= 1;//TODO add modifications                    
                }
                /*
                foreach (Trinity_Gpu.ProteinPrecursorMatch match in pg.Execute(proteinMasses, 0.00005, 10000000))//TODO compute correct tolerance window
                {
                    int size = match.proteinEndPos - match.proteinStartPos;
                    string str = sequence.Substring(match.proteinStartPos, size);
                    if (size >= dbOptions.MinimumPeptideLength)
                    {
                        sequencesTmp[str] -= 1;//TODO add modifications
                    }
                }//*/

                foreach (int val in sequencesTmp.Values)
                    if (val != 0)
                        return false;//*/
            }
            //pg.Dispose();
            return true;
        }
    }
}
