/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.IO;
using System.Threading;
using System.Threading.Tasks;

namespace Trinity
{

    /// <summary>
    /// Regroups method used to map digested peptide sequences to one or more spectrum
    /// </summary>
    public class DBSearcher
    {
        public DBOptions options;
        public DBSearcher(DBOptions options)
        {
            this.options = options;
        }

        //public ConcurrentDictionary<string, List<Protein>> DicOfProteins = new ConcurrentDictionary<string, List<Protein>>();         
        private void ComputePSMs(Query query, Peptide modified_peptide)
        {
            PeptideSpectrumMatch psm = new PeptideSpectrumMatch(query, modified_peptide, options);

            if (psm.MatchingProducts > 1)
            {
                //DicOfProteins.GetOrAdd(psm.Peptide.BaseSequence, new List<Protein>()).Add(psm.Peptide.Parent);
//                if (!DicOfProteins.ContainsKey(psm.Peptide.BaseSequence))
//                    DicOfProteins.AddOrUpdate(psm.Peptide.BaseSequence, new List<Protein>());
//                DicOfProteins[psm.Peptide.BaseSequence].Add(psm.Peptide.Parent);

                //List<PeptideSpectrumMatch> psmsOfQuery = query.psms;

                if (query.psms.Count < options.NbPSMToKeep)
                    query.psms.Add(psm);
                else if (query.psms[options.NbPSMToKeep-1].ProbabilityScore() < psm.ProbabilityScore())
                {
                    for (int i = 0; i < query.psms.Count; i++)
                        if (query.psms[i].ProbabilityScore() <= psm.ProbabilityScore())
                        {
                            query.psms.Insert(i, psm);
                            break;
                        }
                    if (query.psms.Count > options.NbPSMToKeep)
                        query.psms.RemoveAt(options.NbPSMToKeep-1);
                }//*/
                //TODO check optimal number of PSM to store
                /*
                //TODO Reactivate this to save on memory space                
                if (psmsOfPrecursor.Count == 0 || psmsOfPrecursor[0].MaxQuantScore() == psm.MaxQuantScore())
                    psmsOfPrecursor.Add(psm);
                else
                    if (psm.MaxQuantScore() > psmsOfPrecursor[0].MaxQuantScore())
                    {
                        psmsOfPrecursor.Clear();
                        psmsOfPrecursor.Add(psm);
                    }
                //*/                
            }
        }

        /// <summary>
        /// Maps all peptide sequences to potential precursors and spectrum
        /// </summary>
        /// <param name="queries"></param>
        /// <param name="fittingPeptides"></param>
        /// <param name="previousProteins"></param>
        /// <returns></returns>
        public Precursors Search(Queries queries, IEnumerable<Tuple<Peptide, int>> fittingPeptides)
        {
            queries.dbOptions.ConSole.WriteLine("Mapping " + queries.Count + " queries to the digested proteome ... ");

            long nbQueryConsidered = 0;
            Parallel.ForEach<Tuple<Peptide, int>>(fittingPeptides, (Tuple<Peptide, int> hit) =>
            //foreach (Tuple<Peptide, int> hit in fittingPeptides)
            {
                int indexPrecursor = hit.Item2;
                double maximumMass = MassTolerance.MzTop(hit.Item1.MonoisotopicMass, options.precursorMassTolerance);
                double minimumMass = MassTolerance.MzFloor(hit.Item1.MonoisotopicMass, options.precursorMassTolerance);
                
                if (indexPrecursor < queries.Count && queries[indexPrecursor].precursor.Mass >= minimumMass)
                {
                    while (indexPrecursor < queries.Count && queries[indexPrecursor].precursor.Mass <= maximumMass)
                    {

                        lock (queries[indexPrecursor].psms)
                        {
                            //if (low_index < Count && this[low_index].precursor.Mass >= minimum_precursor_mass)
                            //  foreach (Query query in queries.GetQueryInMassRange(modified_peptide.MonoisotopicMass, options.precursorMassTolerance))
                            //{
                            //Target (or decoy with enzyme digests)
                            ComputePSMs(queries[indexPrecursor], hit.Item1);

                            //Decoy if NoEnzyme digest
                            if (options.DecoyFusion)
                                ComputePSMs(queries[indexPrecursor], hit.Item1.Reverse());

                            indexPrecursor++;

                            //foreach (Precursor isotope in query.precursor.Isotopes)                        
                            //    ComputePSMs(query, modified_peptide, isotope.MassShift, previousProteins);
                        }
                    }

                    nbQueryConsidered += indexPrecursor - hit.Item2;
                }
                else
                    options.ConSole.WriteLine("WTF####");
            });

            //Push PSMs to Precursor
            foreach (Query query in queries)
            {
                query.precursor.psms_AllPossibilities.AddRange(query.psms);//No MERGE
                /*
                //Push PSMs to precursors
                if (query.precursor.psms_AllPossibilities.Count == 0)
                    query.precursor.psms_AllPossibilities.AddRange(query.psms);
                else
                {
                    //Merge common entries
                    foreach (PeptideSpectrumMatch psmQuery in query.psms)
                    {
                        bool isNew = true;
                        foreach (PeptideSpectrumMatch psmPrecursor in query.precursor.psms_AllPossibilities)
                        {
                            if (psmPrecursor != psmQuery && psmPrecursor.Peptide == psmQuery.Peptide)//Peptide object should be the same 
                            {
                                psmPrecursor.Merge(psmQuery, options);
                                isNew = false;
                                break;
                            }
                        }
                        if (isNew)
                            query.precursor.psms_AllPossibilities.Add(psmQuery);
                    }
                }//*/
            }
            //PeptideSpectrumMatches allPsms = new PeptideSpectrumMatches();
            int nbAssignedPrecursor = 0;
            foreach (Precursor precursor in queries.Precursors)
                if (precursor.psms_AllPossibilities.Count > 0)
                {
                    nbAssignedPrecursor++;
                //    allPsms.AddRange(precursor.psms);
                }

            //TODO Check the impact of this approach!!!!
            //List<PeptideSpectrumMatch> okPSMs = allPsms.ComputeAtFDR(0.05);//TODO Add parameter for this value            
            //Dictionary<PeptideSpectrumMatch, int> dicOfPsm = new Dictionary<PeptideSpectrumMatch, int>();
            //foreach(PeptideSpectrumMatch match in okPSMs)
            //    dicOfPsm.Add(match, 0);
            //foreach (Precursor precursor in spectra.Precursors)
            //{
            //    List<PeptideSpectrumMatch> newList = new List<PeptideSpectrumMatch>();
            //    foreach (PeptideSpectrumMatch psm in precursor.psms)
            //        if (dicOfPsm.ContainsKey(psm))
            //            newList.Add(psm);
            //    precursor.psms = newList;
            //}

            int nbAssignedQuery = 0;
            foreach (Query query in queries)
                if(query.psms.Count > 0)
                    nbAssignedQuery++;
            options.ConSole.WriteLine(nbAssignedQuery + " queries matched [" + nbAssignedPrecursor + " precursors] out of " + nbQueryConsidered + " psm computed");
                        
            return queries.Precursors;
        }
    }
}