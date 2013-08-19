/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Trinity.MaxQuant
{
    public interface IRawFile
    {
        int MS1Count { get; }
        int NumberOfMS1MassRanges { get; }
        int MS2Count { get; }
        double[] GetMS1MassRange(int i);
        byte GetMS1MassRangeIndex(int i);
        double[] GetMS1TimeSpan(int i);
        Spectrum GetMS1Spectrum(int j, bool subtractBackground, int quantile);
        /*Spectrum GetMS2Spectrum(int i);
        int GetScanNumberFromMs2Index(int i);
        double GetMS2MonoisotopicMz(int i);
        SignalType GetMS2SignalType(int i);//*/
    }
    public class RawFileWrapper: IRawFile
    {
        private Spectra Spectra;
        public RawFileWrapper(Spectra spectra)
        {
            this.Spectra = spectra;
        }

        public int MS1Count { get { return Spectra.MS1s.Count; } }
        public int NumberOfMS1MassRanges { get{ return 1; } }
        public int MS2Count { get { return Spectra.Count; } }
        public double[] GetMS1MassRange(int i)
        {
            return new double[]{Spectra.MS1s[i].MinMz, Spectra.MS1s[i].MaxMz};
        }
        public double[] GetMS1TimeSpan(int i)
        {
            return new double[] { Spectra.MS1s[i].RetentionTimeInMin, Spectra.MS1s[i].RetentionTimeInMin + Spectra.MS1s[i].ScanDuration }
                ;
        }

        public byte GetMS1MassRangeIndex(int i) { return 0; }
        public Spectrum GetMS1Spectrum(int i, bool subtractBackground, int quantile)
        {
            return Spectra.MS1s[i].Peaks;
        }
        /*public Spectrum GetMS2Spectrum(int i);
        public int GetScanNumberFromMs2Index(int i);
        public double GetMS2MonoisotopicMz(int i);
        public SignalType GetMS2SignalType(int i);//*/
    }
}
