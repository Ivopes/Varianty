using DetekceVariant;

var v = new VariantDetector("../../../files/bio.sam", "../../../files/chr17.fa");

var res = v.Calculate();

foreach (var item in res)
{
    Console.WriteLine($"{item.Item1} {item.Item2}");
}