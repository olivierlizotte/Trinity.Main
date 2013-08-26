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

namespace Trinity
{
    /// <summary>
    /// Propheus is the main object of Trinity. At its most basic, it is a modification of the Morpheus 
    /// peptide identification scheme, tailored for no enzyme searches
    /// </summary>
    public class Propheus
    {
        public DBOptions dbOptions;
        public Samples Project;

        //Preloaded lists
        public List<Protein> AllProteins;
        public Dictionary<Sample, Spectra> AllSpectras;
        public Queries AllQueries;
        
        public Propheus(DBOptions dbOptions, Samples Project)
        {
            this.dbOptions = dbOptions;
            this.Project = Project; 
        }

        /// <summary>
        /// Builds the list of proteins to digest, the spectrum to match them to
        /// Also creates the list of queries (some spectrum are used more than once, when multiple 
        /// precursors are found in the mass range of the fragmentation window
        /// </summary>
        public void PrepareForSearch()
        {
            AllProteins             = ReadProteomeFromFasta(dbOptions.FastaDatabaseFilepath, !dbOptions.DecoyFusion);
            AllSpectras             = LoadSpectras();

            AllQueries              = CreateQueries(AllSpectras);
        }

        /// <summary>
        /// Load a Propheus searh object from a Result object
        /// Can be used to load a save state
        /// </summary>
        /// <param name="tmp"></param>
        public void Load(Result tmp)
        {
            AllProteins = ReadProteomeFromFasta(dbOptions.FastaDatabaseFilepath, !dbOptions.DecoyFusion);
            AllQueries = tmp.queries;
            
            AllSpectras = new Dictionary<Sample, Spectra>();
            foreach(Query query in tmp.queries)
            {
                if(!AllSpectras.ContainsKey(query.sample))
                    AllSpectras.Add(query.sample, null);
            }
        }

        /// <summary>
        /// Reads the proteins from a fasta file format
        /// </summary>
        /// <param name="fileName"></param>
        /// <param name="onTheFlyDecoys"> Adds reverse sequnce proteins </param>
        /// <returns></returns>
        public static List<Protein> ReadProteomeFromFasta(string fileName, bool addReverseProteins)
        {
            Console.WriteLine("Reading FASTA file " + fileName + " ... ");
            //Extract Proteins from Fasta file
            List<Protein>  AllProteins = new List<Protein>();
            FileStream protein_fasta_database = new FileStream(fileName, FileMode.Open, FileAccess.Read, FileShare.Read);
            foreach (Protein protein in ProteinFastaReader.ReadProteins(protein_fasta_database, addReverseProteins))
            {
                AllProteins.Add(protein);
            }
            //AllProteins.Sort(Protein.TargetDecoyComparison);
            protein_fasta_database.Close();
            Console.WriteLine("Proteins in fasta file : " + AllProteins.Count / 2 + " [" + AllProteins.Count + "]");
            return AllProteins;
        }        
        
        /// <summary>
        /// Loads spectra from Raw files
        /// </summary>
        /// <returns></returns>
        public Dictionary<Sample, Spectra> LoadSpectras()
        {
            //TODO test compatibility with QExactive, mzML ... other known formats
            AllSpectras = new Dictionary<Sample, Spectra>();
            for (int i = 0; i < Project.Count; i++)
            {
                Sample sample = Project[i];
                string trackFile = dbOptions.OutputFolder + vsCSV.GetFileName_NoExtension(sample.sSDF) + "_Tracks.csv";
                string msmsIonFile = dbOptions.OutputFolder + vsCSV.GetFileName_NoExtension(sample.sSDF) + "_MSMSIons.csv";
                if(dbOptions.LoadSpectraIfFound && System.IO.File.Exists(trackFile)
                                                && System.IO.File.Exists(msmsIonFile))
                {
                    Console.WriteLine("Loading Sectra from " + trackFile + " AND " + msmsIonFile);                    
                    AllSpectras.Add(sample, Spectra.Import(msmsIonFile, trackFile));
                }
                else
                {
                    Console.WriteLine("Loading Sectra " + sample.sSDF);

                    pwiz.CLI.msdata.MSDataFile msFile = new pwiz.CLI.msdata.MSDataFile(sample.sSDF);
                    Spectra spectra = Spectra.Load(msFile, dbOptions, sample.sSDF);
                    spectra.Sort(ProductSpectrum.AscendingPrecursorMassComparison);

                    Console.WriteLine(sample.sSDF + " [" + spectra.Count + " msms scans]");
                    if (dbOptions.SaveMS1Peaks)
                        spectra.ExportTracks(trackFile);

                    if (dbOptions.SaveMSMSPeaks)
                        spectra.ExportMSMS(msmsIonFile);

                    AllSpectras.Add(sample, spectra);
                }           
            }
            return AllSpectras;            
        }
        
        /// <summary>
        /// Creates the list of queries
        /// </summary>
        /// <param name="spectras"></param>
        /// <returns></returns>
        public Queries CreateQueries(Dictionary<Sample, Spectra> spectras)
        {
            Queries AllQueries = new Queries(dbOptions);
            foreach(Sample entry in spectras.Keys)
                AllQueries.GenerateQueries(entry, spectras[entry], spectras[entry].tracks);//TODO Get a list of Tracks, built in LoadSpectra
                        
            return AllQueries;
        }
        
        /// <summary>
        /// Latest version of the search routine. Associates spectrum to digested peptide sequences, 
        /// aligns precursors and fragments, clusters common precursors accross samples, create the list of detected
        /// peptide and protein sequences.
        /// TODO Align maps together (retention time)
        /// </summary>
        /// <param name="queries"></param>
        /// <returns></returns>
        public Result SearchLatestVersion(Queries queries)
        {
            Result result = new Result();
            result.queries = queries;
            result.dbOptions = dbOptions;

            DBSearcher dbSearcher = new DBSearcher(dbOptions);
            Digestion ps = new Digestion(dbOptions);

            if (dbOptions.NoEnzymeSearch)
                result.SetPrecursors(dbSearcher.Search(queries, ps.DigestProteomeOnTheFlyNoEnzyme(AllProteins, queries)));
            else
                result.SetPrecursors(dbSearcher.Search(queries, ps.DigestProteomeOnTheFly(AllProteins, false, queries)));
            Console.WriteLine(result.precursors.Count + " precursors matched !");
            
            PeptideSpectrumMatches allPSMs = new PeptideSpectrumMatches();
            foreach (Precursor precursor in result.precursors)
                foreach (PeptideSpectrumMatch psm in precursor.psms_AllPossibilities)
                    allPSMs.Add(psm);

            //Add all psm possibilities to psms list
            foreach (PeptideSpectrumMatch psm in allPSMs)
                psm.Query.precursor.psms.Add(psm);

            long nbTargets = result.SetPrecursors(result.precursors);

            //Console.WriteLine("Targets before Optimizing Score Ratios : " + nbTargets + " [" + result.matchedPrecursors.Count + "]");

            //allPSMs.OptimizePSMScoreRatios(dbOptions, dbOptions.PSMFalseDiscoveryRate);
            //result.matchedPrecursors.OptimizePSMScoreRatios(dbOptions, dbOptions.PSMFalseDiscoveryRate);
            nbTargets = result.SetPrecursors(result.precursors);
            Console.WriteLine("Targets : " + nbTargets + " [" + result.matchedPrecursors.Count + "]");
            //*/
            //TODO Improve alignment results
            
            Align.AlignPrecursorsByDiff(result, allPSMs);
            nbTargets = result.SetPrecursors(result.precursors);
            Console.WriteLine("Targets after precursor alignment : " + nbTargets + " [" + result.matchedPrecursors.Count + "]");

            Align.AlignProductsByDiff(result, allPSMs);
            nbTargets = result.SetPrecursors(result.precursors);
            Console.WriteLine("Targets after fragment alignment : " + nbTargets + " [" + result.matchedPrecursors.Count + "]");
            //*/
            
            Align.CropPrecursors(result, allPSMs);
            nbTargets = result.SetPrecursors(result.precursors);
            Console.WriteLine("Targets after croping precursors : " + nbTargets + " [" + result.matchedPrecursors.Count + "]");

            Align.CropProducts(result, allPSMs);
            nbTargets = result.SetPrecursors(result.precursors);
            Console.WriteLine("Targets after croping fragments : " + nbTargets + " [" + result.matchedPrecursors.Count + "]");
            //*/
            allPSMs = null;
            long bestTargets = nbTargets;            
            MSSearcher msSearcher = new MSSearcher(dbOptions);
            msSearcher.CumulPsm(result.matchedPrecursors);//TODO Check if its still needed

            //Step 1 : Cluster psms together based on precursor feature //TODO Implement ProteoProfile Scoring based clustering
            //Group in clusters
            result.clusters = msSearcher.Search(result.matchedPrecursors);
            //Todo Align retention times 
            //Todo redo clusterization, based on retention time aligned maps

            //Step 2 : Regroup based on peptide sequence (Morpheus code)
            PeptideSearcher pepSearcher = new PeptideSearcher(dbOptions);
            result.peptides = pepSearcher.Search(result.clusters, result.matchedPrecursors, true);
            result.peptideSequences = pepSearcher.Search(result.clusters, result.matchedPrecursors, false);

            //Step 3 : Regroup based on protein sequences (Morpheus code)
            ProteinSearcher protSearcher = new ProteinSearcher(dbOptions, dbSearcher.DicOfProteins);
            result.proteins = protSearcher.SearchLatest(result.peptides, AllProteins);

            UpdatePsmScores(result.proteins);
            nbTargets = result.SetPrecursors(result.precursors);

            Console.WriteLine("Targets before Optimizing Score Ratios : " + nbTargets + " [" + result.matchedPrecursors.Count + "]");
            result.matchedPrecursors.OptimizePSMScoreRatios(dbOptions, dbOptions.PSMFalseDiscoveryRate);            
            nbTargets = result.SetPrecursors(result.precursors);
            Console.WriteLine("Targets after ReOptimizing PSM Score Ratios : " + nbTargets + " [" + result.matchedPrecursors.Count + "]");
            
            //Step 5 : Compute the new number of Targets
            nbTargets = result.SetPrecursors(result.precursors);
            if (nbTargets < bestTargets)
                Console.WriteLine("FAILED to improve PSMs while adding protein and peptide information");
            
            Console.WriteLine(result.matchedPrecursors.Count + " precursors remaining after ProPheus Search!");                        
            return result;
        }//*/

        /// <summary>
        /// Updates scores for the PeptideSpectrumMatches based on newly computed elements (such as protein sequences, peptides, newly computed tolerances)
        /// </summary>
        /// <param name="protein_groups"></param>
        public static void UpdatePsmScores(List<ProteinGroupMatch> protein_groups)
        {
            //Push ProteinScores AND Peptide Score down to PSMs
            protein_groups.Sort(ProteinGroupMatch.AscendingProbabilityScore);
            foreach (ProteinGroupMatch group in protein_groups)
            {
                double proteinScore = group.ProbabilityScore();
                foreach (PeptideMatch peptide in group.PeptideMatches)
                {
                    double peptideScore = peptide.ProbabilityScore();
                    foreach (Cluster cluster in peptide.clusters)
                        foreach (clCondition condition in cluster.conditions)
                            if (condition != null)
                                foreach (clReplicate replicate in condition.replicates)
                                    if (replicate != null)
                                        foreach (Precursor precursor in replicate.precursors)
                                        {
                                            foreach (PeptideSpectrumMatch psm in precursor.psms)
                                            {
                                                if (peptide.peptide.IsSamePeptide(psm.Peptide, true))
                                                {
                                                    if (psm.ProteinScore < proteinScore)
                                                        psm.ProteinScore = proteinScore;

                                                    if (psm.PeptideScore < peptideScore)
                                                        psm.PeptideScore = peptideScore;
                                                }
                                            }
                                            foreach (Precursor isotope in precursor.Isotopes)
                                                foreach (PeptideSpectrumMatch psm in isotope.psms)
                                                {
                                                    if (psm.ProteinScore < proteinScore)
                                                        psm.ProteinScore = proteinScore;

                                                    if (psm.PeptideScore < peptideScore)
                                                        psm.PeptideScore = peptideScore;
                                                }

                                            precursor.psms.Sort(PeptideSpectrumMatch.DescendingOptimizedScoreComparison);
                                            if (precursor.psms.Count == 0)
                                                Console.WriteLine("Precursor with no match");// = precursor;
                                        }
                }
            }
        }
    }
}
