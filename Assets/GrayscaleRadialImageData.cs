using UnityEngine;

public class GrayscaleRadialImageData : IRadialImageData
{
    BesselSystem _besselSystem;
    BesselSeriesData _seriesData; // one channel

    public GrayscaleRadialImageData(BesselSeriesData seriesData, BesselSystem system)
    {
        _seriesData = seriesData;
        _besselSystem = system;
    }

    public Color GetColor(float r01, float phi)
    {
        float res = 0;
        for (int harmonic = 0; harmonic < _seriesData.HarmonicCount; harmonic++)
        {
            float cos = Mathf.Cos(harmonic * phi);
            float sin = Mathf.Sin(harmonic * phi);
            for (int rootIndex = 0; rootIndex < _seriesData.RootsCount; rootIndex++)
            {
                float bessel = _besselSystem.Bessel(harmonic, rootIndex, r01);
                float member = _seriesData.CosFourierCoefficients[harmonic, rootIndex] * cos * bessel;
                member += _seriesData.SinFourierCoefficients[harmonic, rootIndex] * sin * bessel;
                res += member;
            }
        }
        res = 0.5f + 0.5f * res;
        return new Color(res, res, res);
    }
}