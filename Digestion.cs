﻿/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Proteomics.Utilities;
using System.Threading.Tasks;

namespace Trinity
{
    /// <summary>
    /// Methods to digest proteins into peptides, either with an enzyme, or with the noEnzyme, GPU extension
    /// </summary>
    public class Digestion
    {
        private DBOptions options;
        public Digestion(DBOptions dbOptions)
        {
            this.options = dbOptions;
        }

        private static double[] GetMasses(string sequence)
        {
            double[] proteinMasses = new double[sequence.Length];
            for (int i = 0; i < sequence.Length; i++)
                proteinMasses[i] = AminoAcidMasses.GetMonoisotopicMass(sequence[i]);
            return proteinMasses;
        }

        private bool IsDone = true;
        private System.Collections.Concurrent.ConcurrentQueue<Tuple<Peptide, int>> CumulMatch;
        private void RunGPUProteinDigest(DBOptions options, List<Protein> Proteins, Queries listOfQueries)
        {
            List<double> precursorMasses = new List<double>(listOfQueries.Count);
            //double[] precursors = new double[listOfQueries.Count];
            for (int i = 0; i < listOfQueries.Count; i++)
                precursorMasses.Add(listOfQueries[i].precursor.Mass - Constants.WATER_MONOISOTOPIC_MASS);
            //precursors[i] = listOfQueries[i].precursor.Mass - Constants.WATER_MONOISOTOPIC_MASS;
            precursorMasses.Sort();

            //for each protein, build matrix of mass
            int nbProteins = 0;
            Trinity_Gpu.ProteinDigest pg = new Trinity_Gpu.ProteinDigest(precursorMasses.ToArray(), options.MaximumPeptideLength, options.MinimumPeptideLength);
            foreach (Protein protein in Proteins)
            {
                double[] proteinMasses = GetMasses(protein.BaseSequence);

                foreach (Trinity_Gpu.ProteinPrecursorMatch match in pg.Execute(proteinMasses, options.precursorMassTolerance.Value, options.MaximumPeptideMass))//TODO compute correct tolerance window
                {
                    CumulMatch.Enqueue(new Tuple<Peptide, int>(new Peptide(protein, match.proteinStartPos, match.proteinEndPos, 0), match.firstQueryIndex));
                    //yield return new Tuple<Peptide, int> new Peptide(protein, match.proteinStartPos, match.proteinEndPos, 0);//TODO add modifications
                }

                nbProteins++;
            }
            pg.Dispose();
        }

        /// <summary>
        /// Runs the NoEnzyme protein digestion on the GPU
        /// </summary>
        /// <param name="Proteins"></param>
        /// <param name="listOfQueries"></param>
        /// <returns></returns>
        public IEnumerable<Tuple<Peptide, int>> DigestProteomeOnTheFlyNoEnzyme(List<Protein> Proteins, Queries listOfQueries)
        {
            IsDone = false;
            CumulMatch = new System.Collections.Concurrent.ConcurrentQueue<Tuple<Peptide, int>>();

            Task.Factory.StartNew(() =>
            {
                RunGPUProteinDigest(options, Proteins, listOfQueries);
                IsDone = true;
            });

            Tuple<Peptide, int> item;
            while (!IsDone)
            {
                if (CumulMatch.TryDequeue(out item))
                    yield return item;
            }

            while (CumulMatch.TryDequeue(out item))
            {
                yield return item;
            }
        }

        public IEnumerable<Tuple<Peptide, int>> DigestProteomeOnTheFly(List<Protein> proteins, bool allowSNP, Queries AllQueries)
        {
            //HashSet<string> TargetPeptides = new HashSet<string>();
            //Dictionary<string, int> TargetPeptides = new Dictionary<string, int>();
            //Digest proteins and store peptides in a Dictionnary
            //Does not fit in memory -> 360 Go ....
            //dicOfPeptideSequences = new Dictionary<string, List<Protein>>();
            //double minimumMonoisotopicPeakOffset = dbOptions.precursorMonoisotopicPeakCorrection ? dbOptions.minimumPrecursorMonoisotopicPeakOffset : 0;
            //double maximumMonoisotopicPeakOffset = dbOptions.precursorMonoisotopicPeakCorrection ? dbOptions.maximumPrecursorMonoisotopicPeakOffset : 0;
            foreach (Peptide peptide in ProteinSearcher.ProteinDigest(options, proteins, allowSNP))
            {
                int firstIndex = AllQueries.BinarySearch(MassTolerance.MzFloor(peptide.MonoisotopicMass, options.precursorMassTolerance));
                if (firstIndex >= 0 && firstIndex < AllQueries.Count)
                    yield return new Tuple<Peptide, int>(peptide, firstIndex);
                //foreach (Peptide peptide in ProteinSearcher.ProteinDigestNoEnzyme(dbOptions, proteins, AllQueries))
                //if (!TargetPeptides.Contains(peptide.BaseSequence))
                //{
                foreach (Peptide modPeptide in peptide.GetVariablyModifiedPeptides(options.variableModifications, options.maximumVariableModificationIsoforms))
                {
                    firstIndex = AllQueries.BinarySearch(MassTolerance.MzFloor(modPeptide.MonoisotopicMass, options.precursorMassTolerance));
                    if (firstIndex >= 0 && firstIndex < AllQueries.Count)
                        yield return new Tuple<Peptide, int>(modPeptide, firstIndex);
                }

                //TODO check if this favors targets over decoys since proteins are sorted target->deco
                //    TargetPeptides.Add(peptide.BaseSequence);
                //}            
            }
        }
    }
}