using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Proteomics.Utilities;

namespace Trinity.UnitTest
{
    public class Assignment
    {
        public class Variable
        {
            public string name;
            public double value;
            public double normValue;
            public DateTime time;
            public Variable(string csvLine)
            {
                string[] splits = csvLine.Split(',');
                time = DateTime.Parse(splits[0]);
                name = splits[1];
                value = double.Parse(splits[2]);
            }
            public Variable(DateTime t, string n, double v)
            {
                time = t;
                name = n;
                value = v;
            }
        }

        public static void ExportAllVariables(string fileName, Dictionary<DateTime, List<Variable>> DicOfTime, Dictionary<string, List<Variable>> DicOfVar, Dictionary<string,double> DicOfCorr)
        {
            vsCSVWriter writer = new vsCSVWriter(fileName);
            string title = "Date";
            string corr = "";
            foreach (string key in DicOfVar.Keys)
            {
                title += "," + key;
                corr += "," + DicOfCorr[key];
            }
            writer.AddLine(title);
            writer.AddLine(corr);

            foreach (DateTime key in DicOfTime.Keys)
            {
                string newLine = key.ToString();
                foreach (string nameToMatch in DicOfVar.Keys)
                {
                    string varVal = ",";
                    foreach (Variable v in DicOfTime[key])
                        if (v.name.CompareTo(nameToMatch) == 0)
                        {
                            if (v.name.CompareTo("price") == 0)
                                varVal = "," + v.value;
                            else
                                varVal = "," + v.normValue;
                        }
                    newLine += varVal;
                }
                writer.AddLine(newLine);
            }
            writer.WriteToFile();
        }

        public static void NormalizeVariables(Dictionary<string, List<Variable>> DicOfVar)
        {
            foreach (string key in DicOfVar.Keys)
            {
                double min = double.MaxValue;
                double max = double.MinValue;
                foreach (Variable v in DicOfVar[key])
                {
                    if (v.value < min)
                        min = v.value;
                    if (v.value > max)
                        max = v.value;
                }
                foreach (Variable v in DicOfVar[key])
                    v.normValue = (v.value - min) / (max - min);
            }
        }

        public static void InterpolateMissingValues(string varName, Dictionary<DateTime, List<Variable>> dicOfTime, Dictionary<string, List<Variable>> dicOfVar)
        {
            List<double> times = new List<double>();
            List<double> values = new List<double>();
            double earliestValue = double.MaxValue;
            foreach(Variable v in dicOfVar[varName])
            {
                times.Add(v.time.Ticks);
                values.Add(v.value);
                if (v.time.Ticks < earliestValue && v.value > 0)
                    earliestValue = v.time.Ticks;
            }
            var method = MathNet.Numerics.Interpolation.Interpolate.Common(times, values);
            foreach (DateTime key in dicOfTime.Keys)
            {
                List<Variable> list = dicOfTime[key];
                bool found = false;
                foreach (Variable v in list)
                    if (v.name.CompareTo(varName) == 0)
                        found = true;
                if (!found)
                {
                    if (key.Ticks <= earliestValue)
                        list.Add(new Variable(key, varName, 0));
                    else
                        list.Add(new Variable(key, varName, method.Interpolate(key.Ticks)));
                }
            }
        }

        public static List<double> GetArrayofNormed(string name, Dictionary<DateTime, List<Variable>> dicOfTime, Dictionary<string, List<Variable>> dicOfVar)
        {
            List<double> newList = new List<double>();
            foreach (List<Variable> list in dicOfTime.Values)
            {
                foreach (Variable variable in list)
                {
                    if (variable.name.CompareTo(name) == 0)
                        newList.Add(variable.normValue);
                }
            }
            return newList;
        }

        public static void Launch(IConSol console)
        {
            vsCSV csv = new vsCSV(@"C:\Users\caronlio\Downloads\Via.Science.Pre.Interview.Assignment.Data.2013.10.18.csv");
            Dictionary<DateTime, List<Variable>> DicOfTime = new Dictionary<DateTime, List<Variable>>();
            Dictionary<string, List<Variable>> DicOfVar = new Dictionary<string, List<Variable>>();

            //Data sorted based on date
            for (int i = 1; i < csv.LINES_LIST.Count; i++)
            {
                Variable tmpVar = new Variable(csv.LINES_LIST[i]);
                if (!DicOfTime.ContainsKey(tmpVar.time))
                    DicOfTime.Add(tmpVar.time, new List<Variable>());
                DicOfTime[tmpVar.time].Add(tmpVar);

                if (!DicOfVar.ContainsKey(tmpVar.name))
                    DicOfVar.Add(tmpVar.name, new List<Variable>());
                DicOfVar[tmpVar.name].Add(tmpVar);
            }
            
            foreach(string name in DicOfVar.Keys)
                InterpolateMissingValues(name, DicOfTime, DicOfVar);

            //Rebuild DicOfVar
            DicOfVar.Clear();
            foreach (List<Variable> list in DicOfTime.Values)
            {
                foreach (Variable variable in list)
                {
                    if (!DicOfVar.ContainsKey(variable.name))
                        DicOfVar.Add(variable.name, new List<Variable>());
                    DicOfVar[variable.name].Add(variable);
                }
            }

            //Compute Normalized values
            NormalizeVariables(DicOfVar);            

            //Foreach variable, compare correlation with the "price" variable
            List<double> prices = GetArrayofNormed("price", DicOfTime, DicOfVar);
            Dictionary<string, double> DicOfCorrelation = new Dictionary<string,double>();
            foreach(string name in DicOfVar.Keys)
            {
                List<double> normedVals = GetArrayofNormed(name, DicOfTime, DicOfVar);
                double corr = MathNet.Numerics.Statistics.Correlation.Pearson(prices, normedVals);
                if(name.CompareTo("price") == 0)
                    Console.WriteLine("test");
                DicOfCorrelation.Add(name, corr);
            }

            //Prediction
            vsCSVWriter output = new vsCSVWriter(@"C:\_IRIC\predictions.csv");
            output.AddLine("Time,Price,Prediction");
            foreach (DateTime time in DicOfTime.Keys)
            {
                double pred = 0;
                foreach (string name in DicOfCorrelation.Keys)
                {
                    if(name.CompareTo("price") != 0)
                        foreach (Variable v in DicOfTime[time])
                            if (v.name.CompareTo(name) == 0)
                                pred += DicOfCorrelation[name] * v.normValue;
                }
                pred *= 100000;
                output.AddLine(pred.ToString());
            }
            output.WriteToFile();

            //Export a csv of the varialbes, ordered by date
            ExportAllVariables(@"C:\_IRIC\assignOut.csv", DicOfTime, DicOfVar, DicOfCorrelation);

            Console.WriteLine("Done!");
        }
    }
}
