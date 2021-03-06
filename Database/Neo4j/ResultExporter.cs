﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RestSharp;


namespace Trinity.Database.Neo4j
{
    public static class ResultsToNeo4j
    {
        public static string Export(Result rez)
        {
            RestClient client = new RestClient("http://localhost:7474");

            RestRequest request = new RestRequest("db/data/node", Method.POST);
            request.AddParameter("Type", "Result");
            request.AddParameter("TimeStamp", DateTime.Now.ToLongDateString());
            //request.AddParameter("id", "node"); // adds to POST or URL querystring based on Method
            IRestResponse response = client.Execute(request);
            string[] splits = response.Content.Split('\r');
            string self = "";
            for(int i = 0; i < splits.Length; i++)
                if (splits[i].StartsWith("  \"self\" :"))
                {
                    self = splits[i].Substring(12, splits[i].Length - 12 - 2);
                    break;
                }



            return self;
            //Console.WriteLine(response.Content);
            //return 5;
        }
    }
            /*
            client.PostAsync(request, null);
            request.Po
            client.PostAsync(new RestRequest(
            var response = await client.SendAsync(new HttpRequestMessage(HttpMethod.Post, "cypher/") { Content = new StringContent("{\"query\":\"create me\"}", Encoding.UTF8, "application/json") });
            if(response.StatusCode != HttpStatusCode.OK)
                Console.WriteLine("Not ok");

    for (int j = 0; j < iterations; j++)
    {
        DateTime start = DateTime.Now;
        for (int i = 0; i < amount; i++)
        {
            var response = await client.SendAsync(new HttpRequestMessage(HttpMethod.Post, "cypher/") { Content = new StringContent("{\"query\":\"create me\"}", Encoding.UTF8, "application/json") });
            if(response.StatusCode != HttpStatusCode.OK)
                Console.WriteLine("Not ok");
        }
        TimeSpan timeTaken = DateTime.Now - start;
        Console.WriteLine("took {0}ms", timeTaken.TotalMilliseconds);
    }
            RestClient restClient = new RestClient("//localhost:7474/db/data/node");
            RestRequest restRequest = new RestRequest();
            //restRequest.
            //restClient.Post(
        }
    }//*/

    /// <summary>
    /// Entity base class - includes all properties that should be common among all entities
    /// </summary>
    /*public class Entity
    {
        /// <summary>
        /// Id automatically generated on construction of an object
        /// </summary>
        public Guid Id { get; set; }

        public Entity()
        {
            Id = Guid.NewGuid();
        }
    }//*/
        
    public class ResultExporter : IDisposable
    {
        public void Dispose()
        {
            clientNode = null;
        }
        //private BindingFlags fieldFlags = BindingFlags.Public | BindingFlags.Instance | BindingFlags.GetField;
        //private BindingFlags propertyFlags = BindingFlags.Public | BindingFlags.Instance | BindingFlags.GetProperty;
        private string strGraphList = typeof(IGraphML_List).FullName;
        private string strGraphNode = typeof(IGraphML_Node).FullName;
        
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
        
        Dictionary<long, IGraphML_Node> dic = new Dictionary<long, IGraphML_Node>();
        //GraphClient client;
        RestClient clientNode;
        Result results;
        public ResultExporter(Result resultToExport)
        {
            //client = new GraphClient(new Uri("http://localhost:7474/db/data"));
            //client.Connect();
            //CreateIndexes(client);
            clientNode= new RestClient("//localhost:7474/db/data/node");
        }

        public long AddNode(IGraphML_Node node)
        {
            long id = 0;
            try
            {
                //var nodeRef = client.Create(new Neo4jNode());
                //client.Update(nodeRef,
                //var nodeRef = client.Create(node.GetType());
                //id = nodeRef.Id;
                Type type = node.GetType();
                            /*
                Dictionary<IGraphML_Node, string> relatedVar = new Dictionary<IGraphML_Node,string>();
                Dictionary<IGraphML_List, string> lists = new Dictionary<IGraphML_List, string>();

                StringBuilder strB = new StringBuilder();
                strB.Append("Type : '" + type.FullName + "'");
                //Dictionary<string, string> properties = new Dictionary<string, string>();
                //properties.Add("Type", type.FullName);

                //Properties
                PropertyInfo[] propertyInfos = type.GetProperties(propertyFlags);
                foreach (PropertyInfo property in propertyInfos)
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
                                        if (!string.IsNullOrEmpty(str))
                                            properties.Add(property.Name, ReplaceBadChars(str));
                                    }
                                }
                                catch (Exception ex)
                                {
                                        dbOptions.ConSole.WriteLine(ex.Message);
                                        dbOptions.ConSole.WriteLine(ex.StackTrace);
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
                                if (!string.IsNullOrEmpty(str))
                                    properties.Add(field.Name, ReplaceBadChars(str));
                            }
                        }
                }

                var request = new RestRequest("application/json", Method.POST);
                //request.AddParameter("text/json",,ParameterType.RequestBody);
                //request.AddBody(properties,

                foreach (IGraphML_Node linkedNode in relatedVar.Keys)
                {
                    long id = AddNode(linkedNode, dic);
                    AddRelation(node.ID, id, relatedVar[linkedNode]);
                }

                foreach (IGraphML_List list in lists.Keys)
                {
                    foreach (GraphML_Node item in ((System.Collections.IEnumerable)list))
                    {
                        long id = AddNode(item, dic);
                        AddRelation(node.ID, id, lists[list]);
                    }
                }//*/
            }
            catch (Exception ex)
            {
                results.dbOptions.ConSole.WriteLine(ex.Message);
                results.dbOptions.ConSole.WriteLine(ex.StackTrace);
            }
            return id;// node.ID;
        }
    }
}
