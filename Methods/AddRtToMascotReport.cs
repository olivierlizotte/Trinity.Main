using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Proteomics.Utilities;

namespace Trinity.Methods
{
    public static class AddRtToMascotReport
    {
        public static void FromFolderWithRetentionTimeCSV(string mascotReportCSVFile, string folder, string outputCSV)
        {
            string[] files = Directory.GetFiles(folder, "*_RetentionTimes.csv");
            List<vsCSV> RTs = new List<vsCSV>();
            foreach (string file in files)
                RTs.Add(new vsCSV(file));

            vsCSVWriter writer = new vsCSVWriter(outputCSV);

            vsCSV mascotReport = new vsCSV(mascotReportCSVFile);
            int indexScanNumber = -1;
            int indexFileName = -1;
            int indexRetentionTime = -1;
            bool isContent = false;
            for (int i = 0; i < mascotReport.LINES_LIST.Count; i++)
            {
                string line = mascotReport.LINES_LIST[i];
                string[] splits = line.Split(vsCSV._Generic_Separator);
                if (line.Contains("Scan Number"))
                    indexScanNumber = vsCSV.GetColumnIndex(splits, "Scan Number");
                if(line.Contains("FileName"))
                    indexFileName = vsCSV.GetColumnIndex(splits, "FileName");
                if(line.Contains("Pep Elution Time"))
                    indexRetentionTime = vsCSV.GetColumnIndex(splits, "Pep Elution Time");
                if(isContent)
                {
                    string[] strScanSplits = splits[indexScanNumber].Split('-');                    
                    int tmpScan = 0;
                    for (int k = 0; k < strScanSplits.Length; k++)
                        tmpScan += int.Parse(strScanSplits[k]);
                    
                    string file = vsCSV.GetFileName_NoExtension(splits[indexFileName]);
                    string rt = "";
                    for(int j = 0; j < files.Length; j++)
                        if (files[j].Contains(file))
                        {
                            rt = RTs[j].LINES_LIST[tmpScan].Split(vsCSV._Generic_Separator)[1];
                            break;
                        }
                    
                    splits[indexRetentionTime] = rt;
                    line = vsCSV.Concatenate(splits, ",");
                }
                if (indexScanNumber >= 0 && indexScanNumber < splits.Length &&
                    indexFileName >= 0 && indexFileName < splits.Length &&
                    indexRetentionTime >= 0 && indexRetentionTime < splits.Length)
                {
                    isContent = true;
                }

                writer.AddLine(line);
            }
            writer.writeToFile();
        }
    }
}
