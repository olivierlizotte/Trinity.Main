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
using System.Reflection;

namespace Trinity
{
    /// <summary>
    /// All objects inheriting this interface will be stored in the save state GraphML file
    /// </summary>
    public interface IGraphML_Node
    {
        ulong ID { get; set; }
        //void WriteNode(StreamWriter sw);
        //void WriteEdges(StreamWriter sw);
    }
    public abstract class GraphML_Node : IGraphML_Node
    {
        private ulong _ID = 0;
        public ulong ID { get { return _ID; } set { _ID = value; } }
    }
    
    public interface IGraphML_List : IGraphML_Node
    {

    }
    
    public class GraphML_List<T> : List<T>, IGraphML_List where T : IGraphML_Node
    {
        private ulong _ID = 0;
        public ulong ID { get { return _ID; } set { _ID = value; } }
        public GraphML_List(T[] list) : base(list) { }
        public GraphML_List(int size) : base(size) { }
        public GraphML_List(IEnumerable<T> list) : base(list) { }
        public GraphML_List() : base() { }
    }

    /// <summary>
    /// This class uses generic code (from Reflection) to save the state of all objects originating from the Result object.
    /// The file mimicks a GraphML format (without actually being a valid one).
    /// Reflection is slow for building objects, so these methods are used during devlopment, or to store specific unit tests
    /// </summary>
    public class GraphML
    {
        private static ulong ID = 0;
        public static ulong NextID{get{GraphML.ID++; return GraphML.ID;}}

        public Dictionary<object, ulong> GraphLessObjects = new Dictionary<object, ulong>();
        
        public static string PutBackBadChars(string input)
        {
            string result = input;
            result = result.Replace("&amp;", "&");
            result = result.Replace("&quot;", "\"");
            result = result.Replace("&#92;", "\\");
            result = result.Replace("&#47;", "/");
            result = result.Replace("&#58;", ":");
            result = result.Replace("&#59;", ";");
            result = result.Replace("&#60;", "<");
            result = result.Replace("&#62;", ">");
            return result;
        }

        public static string ReplaceBadChars(string input, bool severe = false)
        {
            string result = "";
            foreach (char c in input.ToCharArray())
            {
                switch (c)
                {
                    case '&': result += "&amp;"; break;
                    case '"': result += "&quot;"; break;
                    case '\\': result += "&#92;"; break;
                    case '\n': if (severe) result += "<br/>"; break;
                    case '\r': if (severe) result += "<br/>"; break;
                    case '/': result += "&#47;"; break;
                    case ':': if (severe) result += "&#58;"; else result += c; break;
                    case ';': if (severe) result += "&#59;"; else result += c; break;
                    case '<': if (severe) result += "&#60;"; else result += c; break;
                    case '>': if (severe) result += "&#62;"; else result += c; break;
                    default: result += c; break;
                }
            }
            return result;
        }

        private object ChangeType(string value, Type newType)
        {
            if (newType.IsEnum)
                return Enum.Parse(newType, value);//Enum.ToObject(newType, value);
            else
                return Convert.ChangeType(value, newType);
        }

        public Result ImportGeneric(string fileName)
        {
            ulong headID = 0;
            Dictionary<ulong, IGraphML_Node> dic = new Dictionary<ulong, IGraphML_Node>();
            string line = "";
            //Result result = null;
            try
            {
                using (StreamReader sr = new StreamReader(fileName))
                {
                    line = sr.ReadLine();

                    IGraphML_Node currentObj = null;
                    Type currentType = null;
                    while (line != null && !line.StartsWith("</Graph>"))
                    {
                        line = line.Trim();
                        if (line.StartsWith("<Graph HeadID="))
                            headID = ulong.Parse(line.Split('"')[1]);
                        if (line.StartsWith("<Node"))
                        {
                            string[] header = line.Split('"');//beginning type empty id bracket
                            string type = header[1];
                            ulong id = ulong.Parse(header[3]);
                            currentType = Type.GetType(type);
                            currentObj = (IGraphML_Node)System.Activator.CreateInstance(currentType);
                            currentObj.ID = id;
                            dic.Add(id, currentObj);
                        }
                        if (line.StartsWith("<Relation"))
                        {
                            string[] splits = line.Split('"');//beginning type empty to empty from bracket
                            string varName = splits[1];
                            object destination = dic[ulong.Parse(splits[3])];
                            object source = dic[ulong.Parse(splits[5])];

                            Type sourceType = source.GetType();
                            Type destinationType = destination.GetType();
                            
                            System.Reflection.FieldInfo fi = sourceType.GetField(varName);
                            if (fi == null)
                            {
                                System.Reflection.PropertyInfo pi = sourceType.GetProperty(varName);
                                if (pi != null)
                                {
                                    if (pi.PropertyType.GetInterface(strGraphList) != null)
                                    {
                                        try
                                        {
                                            var add = pi.PropertyType.GetMethod("Add");
                                            add.Invoke(pi.GetValue(source, null), new[] { destination });
                                        }
                                        catch (Exception ex)
                                        {
                                            Console.WriteLine(ex.Message);
                                        }
                                    }
                                    else
                                        pi.SetValue(source, destination, null);
                                }
                            }
                            else
                            {
                                if (fi.FieldType.GetInterface(strGraphList) != null)
                                {
                                    try
                                    {
                                        var add = fi.FieldType.GetMethod("Add");
                                        add.Invoke(fi.GetValue(source), new[] { destination });
                                    }
                                    catch (Exception ex)
                                    {
                                        Console.WriteLine(ex.Message);
                                    }
                                }
                                else
                                    fi.SetValue(source, destination);
                            }
                        }

                        if (!line.StartsWith("<") && currentObj != null)
                        {
                            string[] splits = line.Split('\'');
                            string varName = splits[1];
                            splits = line.Split('"');
                            string value = PutBackBadChars(splits[1]);
                            System.Reflection.FieldInfo fi = currentType.GetField(varName);
                            if (fi == null)
                            {
                                System.Reflection.PropertyInfo pi = currentType.GetProperty(varName);
                                if (pi != null)
                                    pi.SetValue(currentObj, ChangeType(value, pi.PropertyType), null);
                            }
                            else
                                fi.SetValue(currentObj, ChangeType(value, fi.FieldType));//TODO check if it works for enum
                        }
                        line = sr.ReadLine();
                    }
                    sr.Close();//*/
                }//
            }
            catch (System.Exception ex)
            {
                Console.WriteLine("Unexpected error in " + ex.TargetSite + " : " + ex.Message);
                Console.WriteLine("Stack Trace = " + ex.StackTrace);
            }
            Result result = dic[headID] as Result;
            return result;
        }

        //TODO Replace string by more reflection to speed variable access
        public void Export(string fileName, Result result)
        {
            Console.WriteLine("Saving result to a GraphML file : " + fileName);
            try
            {
                using (StreamWriter sw = new StreamWriter(fileName))
                {
                    Dictionary<IGraphML_Node, ulong> dic = new Dictionary<IGraphML_Node, ulong>();
                    result.ID = NextID;

                    sw.WriteLine("<Graph HeadID=\"" + result.ID + "\">");
                    WriteGenericNode(sw, result, dic);
                    sw.WriteLine("</Graph>");
                    sw.Close();
                    /*
                //----------------NODES----------------------------------//

                    sw.WriteLine("  <Node Type=\"" + result.GetType().FullName + "\" ID=\"" + result.ID + "\">");
                    sw.WriteLine("  </Node>");

                    WriteNode(sw, result.precursors, new string[] {  });
                    foreach (Precursor precursor in result.precursors)
                        Export_Precursor(sw, precursor, result.precursors);
                    WriteRelation(sw, result.ID, result.precursors.ID, "precursors");
                    //GraphML.WriteNode(sw, this, new string[] {"Score", "Decoy", "MassShift", "Charge", "Mass"});

                    //Link to Track
                    //Link to PSMs
                    //Link to Isotpoes
                    //Link to OtherCharges
                    //Link to Sample
                

                //--------------RELATIONS--------------------------------//
                    //WriteRelation(sw, result, result.precursors, "Results");
                    //WriteRelations(sw, result.precursors, result.precursors, "Stores");

                    sw.WriteLine("</Graph>");
                    sw.Close();//*/
                }//
            }
            catch (System.Exception ex)
            {
                Console.WriteLine("Unexpected error in " + ex.TargetSite + " : " + ex.Message);
                Console.WriteLine("Stack Trace = " + ex.StackTrace);
            }            
        }

        public BindingFlags fieldFlags = BindingFlags.Public | BindingFlags.Instance | BindingFlags.GetField;
        public BindingFlags propertyFlags = BindingFlags.Public | BindingFlags.Instance | BindingFlags.GetProperty;
        public string strGraphList = typeof(IGraphML_List).FullName;
        public string strGraphNode = typeof(IGraphML_Node).FullName;

        public ulong WriteGenericNode(StreamWriter sw, IGraphML_Node node, Dictionary<IGraphML_Node, ulong> dic)
        {
            if (!dic.ContainsKey(node))
            {
                try
                {
                    if(node.ID <= 0)
                        node.ID = NextID;

                    dic.Add(node, node.ID);

                    Type type = node.GetType();        
                    sw.WriteLine("  <Node Type=\"" + type.FullName + "\" ID=\"" + node.ID + "\">");
                            
                    Dictionary<IGraphML_Node, string> relatedVar = new Dictionary<IGraphML_Node,string>();
                    Dictionary<IGraphML_List, string> lists = new Dictionary<IGraphML_List, string>();

                    //Properties
                    PropertyInfo[] properties = type.GetProperties(propertyFlags);
                    foreach (PropertyInfo property in properties)
                    {
                        if (!"ID".Equals(property.Name) && property.CanWrite && property.GetSetMethod(true).IsPublic)
                        {
                            if (property.PropertyType.GetInterface(strGraphList) != null)
                                lists.Add((IGraphML_List)property.GetValue(node, null), property.Name);
                            else
                                if (property.PropertyType.GetInterface(strGraphNode) != null)
                                    relatedVar.Add((IGraphML_Node)property.GetValue(node, null), property.Name);
                                else
                                {
                                    try
                                    {
                                        object obj = property.GetValue(node, null);
                                        if (obj != null)
                                        {
                                            string str = obj.ToString();
                                            if(!string.IsNullOrEmpty(str))
                                                sw.WriteLine("      '" + property.Name + "' = \"" + ReplaceBadChars(str) + "\"");
                                        }
                                    }
                                    catch (Exception ex)
                                    {
                                        Console.WriteLine(ex.Message);
                                        Console.WriteLine(ex.StackTrace);
                                    }
                                }
                        }
                    }

                    //Fields
                    FieldInfo[] fields = type.GetFields(fieldFlags);
                    foreach (FieldInfo field in fields)
                    {
                        if (field.FieldType.GetInterface(strGraphList) != null)
                            lists.Add((IGraphML_List)field.GetValue(node), field.Name);
                        else
                            if (field.FieldType.GetInterface(strGraphNode) != null)
                                relatedVar.Add((IGraphML_Node)field.GetValue(node), field.Name);
                            else
                            {
                                object obj = field.GetValue(node);
                                if (obj != null)
                                {
                                    string str = obj.ToString();
                                    if(!string.IsNullOrEmpty(str))
                                        sw.WriteLine("      '" + field.Name + "' = \"" + ReplaceBadChars(str) + "\"");
                                }
                            }
                    }
                    sw.WriteLine("  </Node>");

                    foreach (IGraphML_Node linkedNode in relatedVar.Keys)
                    {
                        ulong id = WriteGenericNode(sw, linkedNode, dic);
                        sw.WriteLine("  <Relation Type=\"" + relatedVar[linkedNode] + "\" TO=\"" + id + "\" FROM=\"" + node.ID + "\">");
                        sw.WriteLine("  </Relation>");
                    }

                    foreach (IGraphML_List list in lists.Keys)
                    {
                        foreach (GraphML_Node item in ((System.Collections.IEnumerable)list))
                        {
                            ulong id = WriteGenericNode(sw, item, dic);
                            sw.WriteLine("  <Relation Type=\"" + lists[list] + "\" TO=\"" + id + "\" FROM=\"" + node.ID + "\">");
                            sw.WriteLine("  </Relation>");
                        }
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine(ex.Message);
                    Console.WriteLine(ex.StackTrace);
                }
            }
            return node.ID;
        }
    }
}
