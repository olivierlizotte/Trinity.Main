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

    public class Query : GraphML_Node, ITargetDecoy
    {
        public Sample                       sample;
        public PeptideSpectrumMatches       psms;
        public ProductSpectrum              spectrum;
        public Precursor                    precursor;
        public int                          spectrumIndex;
        public double                       TrapDistance;

        public double Score
        {            get { return ScoreFct(); }        }

        public bool Decoy
        {            get { return precursor.Decoy; }        }

        public bool Target
        {
            get { return !Decoy; }
        }

        public Query()
        {
            this.psms = new PeptideSpectrumMatches();
        }
        public Query(Sample entry, ProductSpectrum spectrum, Precursor precursor, int spectraIndex = -1)
        {
            this.sample     = entry;
            this.psms       = new PeptideSpectrumMatches();
            this.spectrum   = spectrum;
            this.precursor  = precursor;
            this.spectrumIndex = spectraIndex;
            if (spectrum != null && precursor != null)
                this.TrapDistance = Math.Abs(spectrum.PrecursorMZ - precursor.Track.MZ);
        }

        public double ScoreFct(Peptide peptide = null)
        {
            return precursor.OptimizedScore(peptide);
        }
        /*
        public double Intensity(Peptide peptide)
        {
            return precursor.GetIntensity(peptide);
        }//*/
    }
}
