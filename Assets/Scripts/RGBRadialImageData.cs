using UnityEngine;

public class RGBRadialImageData : IRadialImageData
{
    BesselSystem _besselSystem;
    BesselSeriesData _redSeriesData;
    BesselSeriesData _greenSeriesData;
    BesselSeriesData _blueSeriesData;

    public RGBRadialImageData(BesselSeriesData redData, BesselSeriesData greenData, BesselSeriesData blueData, BesselSystem system)
    {
        _redSeriesData = redData;
        _greenSeriesData = greenData;
        _blueSeriesData = blueData;
        _besselSystem = system;
    }

    public Color GetColor(float r01, float phi)
    {
        float red = _redSeriesData.Eval(r01, phi, _besselSystem);
        red = 0.5f + 0.5f * red;
        float green = _greenSeriesData.Eval(r01, phi, _besselSystem);
        green = 0.5f + 0.5f * green;
        float blue = _blueSeriesData.Eval(r01, phi, _besselSystem);
        blue = 0.5f + 0.5f * blue;
        return new Color(red, green, blue);
    }
}