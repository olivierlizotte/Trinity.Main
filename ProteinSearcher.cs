/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using Proteomics.Utilities;

namespace Trinity
{
    /// <summary>
    /// list of Groups of similar protein sequences
    /// </summary>
    public class ProteinGroupMatches : GraphML_List<ProteinGroupMatch>
    {
        public ProteinGroupMatches(IEnumerable<ProteinGroupMatch> list) : base(list) { }
        public ProteinGroupMatches() { }
        private FDRizer<ProteinGroupMatch> uptimizer;
        public List<ProteinGroupMatch> ComputeAtFDR(double desired_fdr, bool displayValues = true)
        {
            if (uptimizer == null)
            {
                List<Comparison<ProteinGroupMatch>> sorts = new List<Comparison<ProteinGroupMatch>>();
                sorts.Add(CompareMorpheusSummedScore);
                sorts.Add(ComparePeptideMatchesCount);
                uptimizer = new FDRizer<ProteinGroupMatch>(this, sorts, null);
            }
            else
                uptimizer.ReStart();
            return uptimizer.Launch(desired_fdr, displayValues);
        }
        public static int CompareMorpheusSummedScore(ProteinGroupMatch left, ProteinGroupMatch right)
        {
            return -left.SummedMorpheusScore.CompareTo(right.SummedMorpheusScore);
        }
        public static int ComparePeptideMatchesCount(ProteinGroupMatch left, ProteinGroupMatch right)
        {
            return -left.PeptideMatches.Count.CompareTo(right.PeptideMatches.Count);
        }
    }

    /// <summary>
    /// Group of similar protein sequences, that this experiment cannot differentiate
    /// </summary>
    public class ProteinGroupMatch : GraphML_Node, ITargetDecoy
    {
        public GraphML_List<Protein> Proteins;
        public ProteinGroupMatch()
            : base()
        {
            PeptideMatches = new GraphML_List<PeptideMatch>();
            Proteins = new GraphML_List<Protein>();
        }

        public bool Target { get { return !Decoy; } }

        public bool Decoy
        {
            get
            {
                int nbDecoy = 0;
                foreach (PeptideMatch match in PeptideMatches)//Is this method adequate? SummedScore comparison?
                    if (match.Decoy)
                        nbDecoy++;
                return nbDecoy > PeptideMatches.Count - nbDecoy;
            }
        }

        public GraphML_List<PeptideMatch> PeptideMatches { get; set; }
        
        public HashSet<string> BaseLeucinePeptideSequences
        {
            get
            {
                HashSet<string> base_leucine_peptide_sequences = new HashSet<string>();

                foreach (PeptideMatch match in PeptideMatches)
                    base_leucine_peptide_sequences.Add(match.peptide.BaseSequence);

                return base_leucine_peptide_sequences;
            }
        }
        //TODO Rename this since it is not a MorpheusScore anymore
        public double SummedMorpheusScore
        {
            get
            {
                double summed_morpheus_score = 0.0;

                foreach (PeptideMatch match in PeptideMatches)  // need option to score based on all PSMs rather than unique peptide PSMs?
                    summed_morpheus_score += match.ScoreFct();

                return summed_morpheus_score;
            }
        }
        
        public static int DescendingSummedMorpheusScoreProteinGroupComparison(ProteinGroupMatch left, ProteinGroupMatch right)
        {
            try
            {
                int comparison = -(left.SummedMorpheusScore.CompareTo(right.SummedMorpheusScore));
                if (comparison != 0)
                {
                    return comparison;
                }
                else
                {
                    return left.Target.CompareTo(right.Target);
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
            return 0;
        }

        public static int AscendingSummedMorpheusScoreProteinGroupComparison(ProteinGroupMatch left, ProteinGroupMatch right)
        {
            return left.SummedMorpheusScore.CompareTo(right.SummedMorpheusScore);
        }

        public static readonly string Header = "Protein Description,Protein Sequence,Protein Length,Number of Proteins in Group,Number of Peptide-Spectrum Matches,Protein Sequence Coverage (%),Summed Morpheus Score,Decoy";

        public override string ToString()
        {
            StringBuilder description = new StringBuilder();
            StringBuilder sequence = new StringBuilder();
            StringBuilder length = new StringBuilder();
            if (Proteins.Count > 0)
            {
                foreach (Protein protein in this.Proteins)
                {
                    description.Append(protein.Description.Replace(',', ' ') + " / ");
                    sequence.Append(protein.Sequence + '/');
                    length.Append(protein.Sequence.Length.ToString() + '/');
                }
                description = description.Remove(description.Length - 3, 3);
                sequence = sequence.Remove(sequence.Length - 1, 1);
                length = length.Remove(length.Length - 1, 1);
            }
            StringBuilder sb = new StringBuilder();

            sb.Append(description.ToString() + ',');
            sb.Append(sequence.ToString() + ',');
            sb.Append(length.ToString() + ',');
            sb.Append(Proteins.Count.ToString() + ',');
            sb.Append(PeptideMatches.Count.ToString() + ',');
            StringBuilder sequence_coverage = new StringBuilder();
            if (Proteins.Count > 0)
            {
                foreach (Protein protein in this.Proteins)
                    sequence_coverage.Append((CalculateSequenceCoverage(PeptideMatches, protein) * 100.0).ToString() + '/');
                sequence_coverage = sequence_coverage.Remove(sequence_coverage.Length - 1, 1);
            }
            sb.Append(sequence_coverage.ToString() + ',');
            sb.Append(SummedMorpheusScore.ToString() + ',');
            sb.Append(Decoy.ToString());

            return sb.ToString();
        }

        public double CalculateSequenceCoverage(List<PeptideMatch> matches, Protein protein)
        {
            HashSet<int> covered_residues = new HashSet<int>();
            foreach (PeptideMatch match in matches)
            {
                for (int r = match.peptide.StartResidueNumber; r <= match.peptide.EndResidueNumber; r++)
                {
                    covered_residues.Add(r);
                }
            }
            return (double)covered_residues.Count / protein.Length;
        }
    }

    /// <summary>
    /// Methods to search for protein discovered, based on the list of peptides and clusters found
    /// TODO aggregate peptides into proteins based on the quantification (profile)
    /// </summary>
    public class ProteinSearcher
    {
        DBOptions options;
        Dictionary<string, List<Protein>> DicOfProteins;

        public ProteinSearcher(DBOptions options, Dictionary<string, List<Protein>> dicOfProteins)
        {
            this.options = options;
            this.DicOfProteins = dicOfProteins;
        }

        public static List<ProteinGroupMatch> RetrieveBestProteins(List<ProteinGroupMatch> proteins)
        {
            List<ProteinGroupMatch> filteredProteins = new List<ProteinGroupMatch>();
            int nbDesiredProteins = (int) (proteins.Count * 0.9);

            proteins.Sort(ProteinGroupMatch.DescendingSummedMorpheusScoreProteinGroupComparison);
            foreach (ProteinGroupMatch protein in proteins)
            {
                if (!protein.Decoy)
                {
                    filteredProteins.Add(protein);
                    nbDesiredProteins--;
                    if (nbDesiredProteins < 0)
                        break;
                }
            }
            return filteredProteins;
        }

        public ProteinGroupMatches SearchLatest(List<PeptideMatch> peptides, List<Protein> AllProteins)// Dictionary<string, List<Protein>> dicOfPeptides)
        {
            Dictionary<string, ProteinGroupMatch> peptide_proteins;
            Console.WriteLine("Creating list of proteins ... ");
            peptide_proteins = new Dictionary<string, ProteinGroupMatch>();
            foreach (PeptideMatch match in peptides)
            {
                if (!peptide_proteins.ContainsKey(match.peptide.BaseSequence))
                    peptide_proteins.Add(match.peptide.BaseSequence, new ProteinGroupMatch());

                peptide_proteins[match.peptide.BaseSequence].PeptideMatches.Add(match);
            }

            ProteinGroupMatches protein_groups = new ProteinGroupMatches(peptide_proteins.Values);//proteins_by_description.Values);
            protein_groups.Sort(ProteinGroupMatch.DescendingSummedMorpheusScoreProteinGroupComparison);

            Console.WriteLine("Found " + protein_groups.Count + " proteins.");
            Console.WriteLine("Merging undistinguishable proteins...");
            // merge indistinguishable proteins
            for (int i = 0; i < protein_groups.Count - 1; i++)
            {
                ProteinGroupMatch protein_group = protein_groups[i];

                int j = i + 1;
                while (j < protein_groups.Count)
                {
                    ProteinGroupMatch lower_protein_group = protein_groups[j];

                    if (lower_protein_group.SummedMorpheusScore < protein_group.SummedMorpheusScore)
                    {
                        break;
                    }

                    if (lower_protein_group.BaseLeucinePeptideSequences.SetEquals(protein_group.BaseLeucinePeptideSequences))
                    {
                        protein_group.Proteins.AddRange(lower_protein_group.Proteins);  // should only ever be one protein in the group to add
                        protein_groups.RemoveAt(j);
                    }
                    else
                    {
                        j++;
                    }
                }

                Console.Write("\r{0}%   ", ((100 * i) / protein_groups.Count));
            }

            Console.Write("\r{0}%   ", 100);

            //TODO ITS SUPER ULTRA LONG. Shjorten this!
            // remove subset and subsumable protein groups
            int k = protein_groups.Count - 1;
            while (k >= 1)
            {
                ProteinGroupMatch protein_group = protein_groups[k];
                HashSet<string> protein_group_peptides = new HashSet<string>(protein_group.BaseLeucinePeptideSequences);

                for (int l = 0; l < k; l++)
                {
                    ProteinGroupMatch higher_protein_group = protein_groups[l];

                    protein_group_peptides.ExceptWith(higher_protein_group.BaseLeucinePeptideSequences);
                    if (protein_group_peptides.Count == 0)
                    {
                        break;
                    }
                }

                if (protein_group_peptides.Count == 0)
                {
                    protein_groups.RemoveAt(k);
                }
                k--;
            }

            protein_groups.Sort(ProteinGroupMatch.DescendingSummedMorpheusScoreProteinGroupComparison);
            Console.WriteLine("Found " + protein_groups.Count + " protein groups.");
            return protein_groups;
        }

        public static void Export(string filename, List<ProteinGroupMatch> proteins)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine(ProteinGroupMatch.Header);
            foreach (ProteinGroupMatch group in proteins)
                writer.AddLine(group.ToString());
            writer.writeToFile();
        }

        public static IEnumerable<Peptide> ProteinDigest(DBOptions options, List<Protein> Proteins, bool allowSNP, bool onTheFlyNoEnzymeDecoys)
        {
            int nbProteins = 0;
            bool tooLong = false;
            foreach (Protein protein in Proteins)
            {
                nbProteins++;
                if (tooLong)
                    break;
                List<int> indices = options.DigestionEnzyme.GetDigestionSiteIndices(protein);
                indices.Insert(0, -1);
                indices.Add(protein.Length - 1);

                for (int missed_cleavages = 0; missed_cleavages <= options.ToleratedMissedCleavages; missed_cleavages++)
                {
                    for (int i = 0; i < indices.Count - missed_cleavages - 1; i++)
                    {
                        if (indices[i + missed_cleavages + 1] + 1 - (indices[i] + 1 + 1) + 1 >= options.MinimumPeptideLength)
                        {
                            if (options.initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || indices[i] + 1 != 0 || protein[0] != 'M')
                            {
                                Peptide newPep = new Peptide(protein, indices[i] + 1, indices[i + missed_cleavages + 1], missed_cleavages);
                                yield return newPep;

                                if (allowSNP)
                                    foreach (Peptide possibleSnp in newPep.GetSNPsdPeptides())
                                        yield return possibleSnp;

                                if (!onTheFlyNoEnzymeDecoys)
                                    yield return new Peptide(protein, indices[i + missed_cleavages + 1], indices[i] + 1, missed_cleavages);
                            }

                            if (options.initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && indices[i] + 1 == 0 && protein[0] == 'M')
                            {
                                Peptide newPep = new Peptide(protein, indices[i] + 1 + 1, indices[i + missed_cleavages + 1], missed_cleavages);
                                yield return newPep;

                                if (allowSNP)
                                    foreach (Peptide possibleSnp in newPep.GetSNPsdPeptides())
                                        yield return possibleSnp;

                                if (!onTheFlyNoEnzymeDecoys)
                                    yield return new Peptide(protein, indices[i + missed_cleavages + 1], indices[i] + 1 + 1, missed_cleavages);
                            }
                        }
                    }
                }
                Console.Write("\r{0}%   ", ((100 * nbProteins) / Proteins.Count));
            }
            Console.Write("\r{0}%   ", 100);
        }      
    }
}