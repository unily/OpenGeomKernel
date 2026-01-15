#include "Geom_CurveOnSurface.hxx"

#include <gp_Pnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <GeomAPI.hxx>
#include <Standard_Type.hxx> // 若需

// 在 .cpp 中定义 RTTI 的实现
IMPLEMENT_STANDARD_RTTIEXT(Geom_CurveOnSurface, Geom_BoundedCurve)

//=======================================================================
//function : 构造
//=======================================================================
Geom_CurveOnSurface::Geom_CurveOnSurface(const Handle(Geom2d_BoundedCurve)& pCurve,
    const Handle(Geom_Surface)& surface) : myPCurve(pCurve)
    , mySurface(surface)
{

}

//=======================================================================
//function : Reverse
//=======================================================================
void Geom_CurveOnSurface::Reverse()
{
    if (!myPCurve.IsNull()) {
        myPCurve = Handle(Geom2d_BoundedCurve)::DownCast(myPCurve->Copy());
        myPCurve->Reverse();
    }
}

//=======================================================================
//function : ReversedParameter
//=======================================================================
Standard_Real Geom_CurveOnSurface::ReversedParameter(const Standard_Real U) const
{
    if (!myPCurve.IsNull())
        return myPCurve->ReversedParameter(U);
    return -U; // 简单默认实现
}

//=======================================================================
//function : FirstParameter
//=======================================================================
Standard_Real Geom_CurveOnSurface::FirstParameter() const
{
    if (!myPCurve.IsNull())
        return myPCurve->FirstParameter();
    return Precision::Infinite(); // or RealFirst()
}

//=======================================================================
//function : LastParameter
//=======================================================================
Standard_Real Geom_CurveOnSurface::LastParameter() const
{
    if (!myPCurve.IsNull())
        return myPCurve->LastParameter();
    return -Precision::Infinite(); // or RealLast()
}

//! Returns the end point of the curve.
gp_Pnt Geom_CurveOnSurface::EndPoint() const {
    const gp_Pnt2d& P2 = myPCurve->EndPoint();
    gp_Pnt P;
    mySurface->D0(P2.X(), P2.Y(), P);
    return P;
}

//! Returns the start point of the curve.
gp_Pnt Geom_CurveOnSurface::StartPoint() const {
    const gp_Pnt2d& P2 = myPCurve->StartPoint();
    gp_Pnt P;
    mySurface->D0(P2.X(), P2.Y(), P);
    return P;
}

//=======================================================================
//function : IsClosed
//=======================================================================
Standard_Boolean Geom_CurveOnSurface::IsClosed() const
{
    if (!myPCurve.IsNull())
        return myPCurve->IsClosed();
    return Standard_False;
}

//=======================================================================
//function : IsPeriodic
//=======================================================================
Standard_Boolean Geom_CurveOnSurface::IsPeriodic() const
{
    if (!myPCurve.IsNull())
        return myPCurve->IsPeriodic();
    return Standard_False;
}

//=======================================================================
//function : Continuity
//=======================================================================
GeomAbs_Shape Geom_CurveOnSurface::Continuity() const
{
    if (!mySurface.IsNull() && !myPCurve.IsNull())
    {
        GeomAbs_Shape cSurf = mySurface->Continuity();
        GeomAbs_Shape cPCrv = myPCurve->Continuity();
        return (cSurf < cPCrv) ? cSurf : cPCrv;
    }
    return GeomAbs_C0;
}

//=======================================================================
//function : IsCN
//=======================================================================
Standard_Boolean Geom_CurveOnSurface::IsCN(const Standard_Integer N) const
{
    if (myPCurve.IsNull())
        return Standard_False;
    return myPCurve->IsCN(N);
}

//=======================================================================
//function : D0
//=======================================================================
void Geom_CurveOnSurface::D0(const Standard_Real U, gp_Pnt& P) const
{
    if (myPCurve.IsNull() || mySurface.IsNull())
    {
        P.SetCoord(0., 0., 0.);
        return;
    }
    gp_Pnt2d uv = myPCurve->Value(U);
    P = mySurface->Value(uv.X(), uv.Y());
}

//=======================================================================
//function : D1
//=======================================================================
void Geom_CurveOnSurface::D1(const Standard_Real U, gp_Pnt& P, gp_Vec& V1) const
{
    if (myPCurve.IsNull() || mySurface.IsNull())
    {
        P.SetCoord(0., 0., 0.);
        V1.SetCoord(0., 0., 0.);
        return;
    }

    gp_Pnt2d uv;
    gp_Vec2d uvD1;
    myPCurve->D1(U, uv, uvD1);

    gp_Pnt p3d;
    gp_Vec dU, dV;
    mySurface->D1(uv.X(), uv.Y(), p3d, dU, dV);

    P = p3d;
    V1 = dU.Multiplied(uvD1.X());
    V1.Add(dV.Multiplied(uvD1.Y()));
}

//=======================================================================
//function : D2
//=======================================================================
void Geom_CurveOnSurface::D2(const Standard_Real U, gp_Pnt& P, gp_Vec& V1, gp_Vec& V2) const
{
    if (myPCurve.IsNull() || mySurface.IsNull())
    {
        P.SetCoord(0., 0., 0.);
        V1.SetCoord(0., 0., 0.);
        V2.SetCoord(0., 0., 0.);
        return;
    }

    gp_Pnt2d uv;
    gp_Vec2d uvD1, uvD2;
    myPCurve->D2(U, uv, uvD1, uvD2);

    gp_Pnt p3d;
    gp_Vec dU, dV, d2Uu, d2Uv, d2Vv;
    mySurface->D2(uv.X(), uv.Y(), p3d, dU, dV, d2Uu, d2Uv, d2Vv);

    // D0
    P = p3d;

    // D1
    V1 = dU.Multiplied(uvD1.X());
    V1.Add(dV.Multiplied(uvD1.Y()));

    // D2
    Standard_Real du = uvD1.X();
    Standard_Real dv = uvD1.Y();
    Standard_Real d2u = uvD2.X();
    Standard_Real d2v = uvD2.Y();

    gp_Vec vPart1 = d2Uu.Multiplied(du * du);
    vPart1.Add(d2Uv.Multiplied(2.0 * du * dv));
    vPart1.Add(d2Vv.Multiplied(dv * dv));

    gp_Vec vPart2 = dU.Multiplied(d2u);
    vPart2.Add(dV.Multiplied(d2v));

    V2 = vPart1;
    V2.Add(vPart2);
}

//=======================================================================
//function : D3
//=======================================================================
void Geom_CurveOnSurface::D3(const Standard_Real U, gp_Pnt& P, gp_Vec& V1, gp_Vec& V2, gp_Vec& V3) const
{
    // 简易: 调用 D2, 把 V3=零
    D2(U, P, V1, V2);
    V3.SetCoord(0., 0., 0.);
    // 若需三阶真实展开，需要 Surface->D3(...) + chain rule
}

//=======================================================================
//function : DN
//=======================================================================
gp_Vec Geom_CurveOnSurface::DN(const Standard_Real U, const Standard_Integer N) const
{

    if (N == 1)
    {
        gp_Pnt p;
        gp_Vec v;
        D1(U, p, v);
        return v;
    }
    else if (N == 2)
    {
        gp_Pnt p;
        gp_Vec v1, v2;
        D2(U, p, v1, v2);
        return v2;
    }
    // 更高阶自行实现
    return gp_Vec(0., 0., 0.);
}

//=======================================================================
//function : Transform
//=======================================================================
void Geom_CurveOnSurface::Transform(const gp_Trsf& T)
{
    if (!mySurface.IsNull())
        mySurface->Transform(T);
    // 真正要保持 (u,v) - 3D 对应，还需更多处理
}

//=======================================================================
//function : TransformedParameter
//=======================================================================
Standard_Real Geom_CurveOnSurface::TransformedParameter(const Standard_Real U, const gp_Trsf& /*T*/) const
{
    return U; // 简单实现: 不改参数
}

//=======================================================================
//function : ParametricTransformation
//=======================================================================
Standard_Real Geom_CurveOnSurface::ParametricTransformation(const gp_Trsf& /*T*/) const
{
    return 1.0; // 简单实现: 返回1
}

//=======================================================================
//function : Copy
//=======================================================================
Handle(Geom_Geometry) Geom_CurveOnSurface::Copy() const
{
    Handle(Geom2d_BoundedCurve) newPCurve;
    if (!myPCurve.IsNull())
        newPCurve = Handle(Geom2d_BoundedCurve)::DownCast(myPCurve->Copy());

    Handle(Geom_Surface) newSurface;
    if (!mySurface.IsNull())
        newSurface = Handle(Geom_Surface)::DownCast(mySurface->Copy());

    Handle(Geom_CurveOnSurface) theCopy = new Geom_CurveOnSurface(newPCurve, newSurface);
    return theCopy;
}

//=======================================================================
//function : DumpJson
//=======================================================================
void Geom_CurveOnSurface::DumpJson(Standard_OStream& theOStream, Standard_Integer /*theDepth*/) const
{
    theOStream << "{ \"Geom_CurveOnSurface\": {\n";
    theOStream << "  \"pCurve\": " << (myPCurve.IsNull() ? "\"NULL\"" : "\"2D curve\"") << ",\n";
    theOStream << "  \"surface\": " << (mySurface.IsNull() ? "\"NULL\"" : "\"3D surface\"") << "\n";
    theOStream << "} }" << std::endl;
}