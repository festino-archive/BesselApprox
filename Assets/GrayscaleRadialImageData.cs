using UnityEngine;

public class GrayscaleRadialImageData : IRadialImageData
{
    BesselSystem _besselSystem;
    BesselSeriesData _graySeriesData; // one channel

    public GrayscaleRadialImageData(BesselSeriesData seriesData, BesselSystem system)
    {
        _graySeriesData = seriesData;
        _besselSystem = system;
    }

    public Color GetColor(float r01, float phi)
    {
        float res = _graySeriesData.Eval(r01, phi, _besselSystem);
        res = 0.5f + 0.5f * res;
        return new Color(res, res, res);
    }
}