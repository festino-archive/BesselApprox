
using UnityEngine;

//Prec is: (Bessel calc settings), number of harmonics, number of Bessel_n roots, integral resolution
public struct RadialImagePrecision
{
    public int HarmonicCount;
    public int RootsCount;

    public RadialImagePrecision(int harmonicCount, int rootsCount)
    {
        HarmonicCount = harmonicCount;
        RootsCount = rootsCount;
    }

    public int IntegralPointsR(int n, int k, int width, int height)
    {
        return Mathf.CeilToInt(Mathf.Sqrt(width * width + height * height));
    }

    public int IntegralPointsPhi(int n, int k, int width, int height, float r01)
    {
        return Mathf.CeilToInt(1 + 2 * Mathf.PI * r01 * Mathf.Sqrt(width * width + height * height));
    }
}

public enum BesselConvertingMode
{
    RGB, GRAYSCALE
}

public struct BesselSeriesData
{
    public int HarmonicCount { get => CosFourierCoefficients.GetLength(0); }
    public int RootsCount { get => CosFourierCoefficients.GetLength(1); }

    public readonly float[,] CosFourierCoefficients;
    public readonly float[,] SinFourierCoefficients;

    public BesselSeriesData(int harmCount, int rootCount)
    {
        CosFourierCoefficients = new float[harmCount, rootCount];
        SinFourierCoefficients = new float[harmCount, rootCount];
    }

    public float Eval(float r01, float phi, BesselSystem _besselSystem)
    {
        float res = 0;
        for (int harmonic = 0; harmonic < HarmonicCount; harmonic++)
        {
            float cos = Mathf.Cos(harmonic * phi);
            float sin = Mathf.Sin(harmonic * phi);
            for (int rootIndex = 0; rootIndex < RootsCount; rootIndex++)
            {
                float bessel = _besselSystem.Bessel(harmonic, rootIndex, r01);
                float member = CosFourierCoefficients[harmonic, rootIndex] * cos * bessel;
                member += SinFourierCoefficients[harmonic, rootIndex] * sin * bessel;
                res += member;
            }
        }
        return res;
    }
}