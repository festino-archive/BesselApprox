using UnityEngine;

public class RGBRadialImageData : IRadialImageData
{
    public Color GetColor(float r, float phi)
    {
        return new Color(r, 1 - r, r);
    }
}