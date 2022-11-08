using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DetekceVariant
{
    internal class VariantDetector
    {
        string[] _seqArr;
        string[] _seqQualityArr;

        int[] _startIndexArr;

        string _refGenom;

        double _h = 0.001;
        double _refErr = 0.01;
        double _base = 1;
        string[] _genotypes = new string[] { "AA", "CC", "GG", "TT", "AC", "AG", "AT", "CG", "CT", "GT" };

        public VariantDetector(string samFile, string refFile)
        {

            ReadSamFile(samFile);

            ReadRefFile(refFile);
        }

        private void ReadRefFile(string refFile)
        {
            double e = GetErrForChar('A');


            _refGenom = File.ReadAllText(refFile);

            _refGenom = _refGenom.Substring(6); // ignore prvni radek
            _refGenom = _refGenom.Replace("\n", "").ToUpper();
        }

        public List<(int, string)> Calculate()
        {
            int[] distinctIndexes = _startIndexArr.Distinct().ToArray();
            List<string[]> variantsResult = new();
            int varIndex = -1;
            List<(int, string)> result = new();

            foreach (int index in distinctIndexes)
            {
                ++varIndex;
                List<string> sequences = new();
                List<string> seqQuality = new();
                for (int i = 0; i < _startIndexArr.Length; i++)
                {
                    if (_startIndexArr[i] == index) 
                    {
                        sequences.Add(_seqArr[i]);
                        seqQuality.Add(_seqQualityArr[i]);
                    }
                }
                var ss = sequences.GroupBy(s => s[0]).ToList();
                int maxLengthOfSequence = sequences.Max(s => s.Length);
                variantsResult.Add(new string[maxLengthOfSequence]);

                for (int j = 0; j < maxLengthOfSequence; j++)
                {
                    List<double> variantProbs = new();
                    foreach (var genotype in _genotypes)
                    {
                        double podminenyResult = 0;
                        for (int k = 0; k < sequences.Count; k++)
                        {
                            string seq = sequences[k];
                            string seqErr = seqQuality[k];
                            if (j >= seq.Length) continue;
                            char c = seq[j];
                            double err = GetErrForChar(seqErr[j]);

                            var d = (double)CalculateProduct(c, genotype[0], genotype[1], err);// * 1000000000;
                            if (podminenyResult * d == 0)
                            {

                            }
                            podminenyResult += Math.Log10(d);

                            if (podminenyResult == 0)
                            {
                                break;
                            }
                        }

                        char refCharr = _refGenom[index - 1 + j];
                        podminenyResult *= (double)GetPforRef(refCharr, genotype, _h, _refErr);
                        variantProbs.Add(podminenyResult);
                    }
                    double first, second;
                    first = second = variantProbs[0];
                    int fI = 0, sI = 0;
                    for (int i = 0; i < variantProbs.Count; i++)
                    {
                        double p = variantProbs[i];
                        if (p > first)
                        {
                            second = first;
                            sI = fI;
                            first = p;
                            fI = i;
                        }
                        else if (p > second)
                        {
                            second = p;
                            sI = i; 
                        }
                    }
                    

                    double lod = first - second;
                    lod = Math.Abs(lod);

                    char refChar = _refGenom[index - 1 + j];
                    if (lod > 0.0001)
                    {
                        string var = GetVariantForRef(refChar, _genotypes[fI]);
                        variantsResult[varIndex][j] = string.Empty;
                        if (!string.IsNullOrEmpty(var))
                        {
                            var name = GetVarName(refChar, var);
                            variantsResult[varIndex][j] = $"{refChar} - {var} ({name})";
                        }
                    }
                    else
                    {
                        variantsResult[varIndex][j] = "";
                        continue;
                        var name1 = GetVarName(refChar, _genotypes[fI]);
                        var name2 = GetVarName(refChar, _genotypes[sI]);

                        variantsResult[varIndex][j] = $"{refChar} - {_genotypes[fI]} ({name1}) / {_genotypes[sI]} ({name2})";
                    }
                    /*
                    if (first == 0 && second == 0)
                    {
                        variantsResult[varIndex][j] = "";
                        continue;
                    }
                    char refChar = _refGenom[index - 1 + j];
                    if (second == 0)
                    {
                        string var = GetVariantForRef(refChar, _genotypes[fI]);
                        variantsResult[varIndex][j] = string.Empty;
                        if (!string.IsNullOrEmpty(var))
                        {
                            var name = GetVarName(refChar, var);
                            variantsResult[varIndex][j] = $"{refChar} - {var} ({name})";
                        }
                    }
                    else
                    {
                        //double lod = Math.Log(first) - Math.Log(second);
                        double lod = first - second;
                        lod = Math.Abs(lod);
                        if (lod > 0.0001)
                        {
                            var name = GetVarName(refChar, _genotypes[fI]);
                            variantsResult[varIndex][j] = $"{refChar} - {_genotypes[fI]} ({name})";
                        }
                        else
                        {
                            variantsResult[varIndex][j] = "";
                            continue;
                            var name1 = GetVarName(refChar, _genotypes[fI]);
                            var name2 = GetVarName(refChar, _genotypes[sI]);

                            variantsResult[varIndex][j] = $"{refChar} - {_genotypes[fI]} ({name1}) / {_genotypes[sI]} ({name2})";
                        }
                    }*/

                    if (!string.IsNullOrEmpty(variantsResult[varIndex][j]))
                    {
                        result.Add((index + j, variantsResult[varIndex][j]));
                    }
                }
            }

            return result;
        }
        private double CalculateProduct(char c, char p1, char p2, double err)
        {
            return CalculateProduct(c, p1, err) / 2 + CalculateProduct(c, p2, err) / 2;
        }
        private double CalculateProduct(char c, char p, double err)
        {
            if (c == p)
            {
                return 1 - err;
            }
            return err / 3;
        }
        private double GetErrForChar(char c)
        {
            double err = Math.Pow(10, - (c - 33) / 10.0);

            return Math.Round(err, 7);
        }
        private string[] GenVariantsHomoRef(char c)
        {
            return new string[] { c.ToString() + c };
        }
        private string[] GenVariantsHomo(char c)
        {
            var res = new List<string>() { "AA", "CC", "GG", "TT" };
            for (int i = res.Count - 1; i >= 0; i--)
            {
                if (res[i].Contains(c)) res.RemoveAt(i);
            }
            return res.ToArray();
        }
        private string[] GenVariantsHeteRef(char c)
        {
            var res = new List<string>() { "A", "C", "G", "T" };
            for (int i = res.Count - 1; i >= 0; i--)
            {
                if (res[i].Contains(c))
                {
                    res.RemoveAt(i);
                    continue;
                }
                res[i] = c + res[i];
            }
            return res.ToArray();
        }
        private string[] GenVariantsHete(char c)
        {
            var res = new List<string>() { "A", "C", "G", "T" };
            for (int i = res.Count - 1; i >= 0; i--)
            {
                if (res[i].Contains(c))
                {
                    res.RemoveAt(i);
                    continue;
                }
            }
            var ress = new string[3];
            ress[0] = res[0] + res[1];
            ress[1] = res[1] + res[2];
            ress[2] = res[0] + res[2];

            return ress.ToArray();
        }
        private double GetPforRef(char c, string genotype, double h, double err)
        {
            if (genotype[0] == c && genotype[1] == c) return 1 - (3 * h) / 2;
            if (genotype[0] == c || genotype[1] == c) return h;
            if (genotype[0] == genotype[1]) return h / 2;
            return h * err;
        }
        private void ReadSamFile(string samFile)
        {
            var lines = File.ReadAllLines(samFile);

            _seqArr = new string[lines.Length - 25]; // zajima nas zacatek na line 26
            _seqQualityArr = new string[_seqArr.Length];
            _startIndexArr = new int[_seqArr.Length];

            for (int i = 25; i < lines.Length; i++)
            {
                var line = lines[i];

                var l = line.Split('\t');

                _seqArr[i - 25] = l[9];
                _seqQualityArr[i - 25] = l[10];
                _startIndexArr[i - 25] = int.Parse(l[3]);

            }

        }
        private string GetVariantForRef(char c, string genotype)
        {
            if (GenVariantsHomoRef(c)[0] == genotype) return "";
            return genotype;
        }
        private string GetVarName(char refC, string genotype)
        {
            if (genotype[0] == refC && genotype[1] == refC) return "Bez varianty";
            if (genotype[0] == refC || genotype[1] == refC) return "Homozygotní";
            if (genotype[0] == genotype[1]) return "Heterozygotní vuci referenci";
            return "Heterozygotní mimo referenci";
        }
    }
}
