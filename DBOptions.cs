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

namespace Trinity
{
    /// <summary>
    /// Stores options related to a sample file or group of files
    /// Values are adjusted during calculations to reflec detected thresholds
    /// </summary>
    public class DBOptions : GraphML_Node
    {
        public int MinimumPrecursorChargeState;
        public int MaximumPrecursorChargeState;
        public int MaximumNumberOfFragmentsPerSpectrum;
        public string FastaDatabaseFilepath;
        public double MaximumPeptideMass;
        public bool DecoyFusion;

        public Protease DigestionEnzyme;
        public int ToleratedMissedCleavages;

        public double PSMFalseDiscoveryRate = 0.05;

        public string OutputFolder;
        public double MinimumPSMScore;
        public int MinimumPeptideLength = 5;
        public int MaximumPeptideLength = 100;
        public Fragments fragments;
        public bool addFragmentLoss = false;
        public bool addFragmentMods = false;
        public bool NoEnzymeSearch = false;
        public bool WriteMaxQuantPeakFile = false;

        public bool SaveMS1Peaks = true;
        public bool SaveMSMSPeaks = true;

        public bool LoadSpectraIfFound = true;

        public double ComputedRetentionTimeDiff = 1.0;//TODO compute this after alignment step, based on common identifications
        public double EffectiveIsolationWindowRatio = 0.9;//TODO optimize this value 
        public double MinimumPrecursorIntensityRatioInIsolationWindow = 0.05;


        public double dProduct = 0.0;//0.01;
        public double dPrecursor = 0.12;//.1;//0.01;
        public double dMatchingProductFraction = 0.45;//0.6;//0.8;//0.26;
        public double dMatchingProduct = 0;//0.6;//0.8;//0.26;
        public double dIntensityFraction = 0.13;//0.8;//0.4;//0.47;
        public double dIntensity = 0;//0.4;//0.47;
        public double dProtein = 0;//0.0002;//0.01;
        public double dPeptideScore = 0.3;//0.008;//0;      

        ///Values kept from the original Morpheus source code
        public InitiatorMethionineBehavior initiatorMethionineBehavior;
        public GraphML_List<Modification> fixedModifications;
        public GraphML_List<Modification> variableModifications;
        public int maximumVariableModificationIsoforms;
        public MassTolerance precursorMassTolerance;
        public MassTolerance productMassTolerance;

        /// <summary>
        /// Parameter less constructor used for save states
        /// </summary>
        public DBOptions()
        {
            fixedModifications = new GraphML_List<Modification>();
            variableModifications = new GraphML_List<Modification>();
            fragments = new Fragments();
        }

        public DBOptions(string fasta)
        {
            //Create with default values
            this.DecoyFusion = true;
            this.FastaDatabaseFilepath = fasta;
            this.MaximumPeptideMass = 10000;
            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            this.DigestionEnzyme = proteases["no enzyme"];// proteases["trypsin (no proline rule)"];
            this.NoEnzymeSearch = true;
            this.ToleratedMissedCleavages = 100;// 3;//determines the length of peptides with no-enzyme option
            this.initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            this.fixedModifications = new GraphML_List<Modification>();
            this.variableModifications = new GraphML_List<Modification>();
            this.maximumVariableModificationIsoforms = 1024;
            
            this.MinimumPrecursorChargeState = 1;
            this.MaximumPrecursorChargeState = 4;
            this.MaximumNumberOfFragmentsPerSpectrum = 400;
            //TODO Add precision to the precursor by reading MS part of file
            this.precursorMassTolerance = new MassTolerance(0.005, MassToleranceUnits.Da);//2.1
            //TODO Add precision to the product masses by reading corresponding MS part of raw file
            this.productMassTolerance = new MassTolerance(0.005, MassToleranceUnits.Da);
            
            this.PSMFalseDiscoveryRate = 0.05;

            this.OutputFolder = @"C:\_IRIC\DATA\Test2";//C:\Documents and Settings\ProteoAdmin\Desktop\AEffacer\Morpheus\Output";
            this.MinimumPSMScore = 0.0001;         
        }
    }
}
