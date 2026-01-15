#ifndef _Geom_CurveOnSurface_HeaderFile
#define _Geom_CurveOnSurface_HeaderFile

#include <Geom_BoundedCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_BoundedCurve.hxx>
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <Standard_Real.hxx>
#include <Standard_Handle.hxx>
#include <Standard_DefineHandle.hxx>
#include <Precision.hxx>

//! 类声明
class Geom_CurveOnSurface : public Geom_BoundedCurve
{
public:
    //! 构造函数
    Standard_EXPORT Geom_CurveOnSurface(const Handle(Geom2d_BoundedCurve)& baseCurve, const Handle(Geom_Surface)& surface);

    //! 翻转方向
    Standard_EXPORT void Reverse() Standard_OVERRIDE;
    //! 计算翻转后参数
    Standard_EXPORT Standard_Real ReversedParameter(const Standard_Real U) const Standard_OVERRIDE;

    //! 返回首参数
    Standard_EXPORT Standard_Real FirstParameter() const Standard_OVERRIDE;
    //! 返回尾参数
    Standard_EXPORT Standard_Real LastParameter() const Standard_OVERRIDE;

    //! Returns the end point of the curve.
    Standard_EXPORT virtual gp_Pnt EndPoint() const Standard_OVERRIDE;

    //! Returns the start point of the curve.
    Standard_EXPORT virtual gp_Pnt StartPoint() const Standard_OVERRIDE;

    //! 是否闭合
    Standard_EXPORT Standard_Boolean IsClosed() const Standard_OVERRIDE;
    //! 是否周期
    Standard_EXPORT Standard_Boolean IsPeriodic() const Standard_OVERRIDE;
    //! 连续性
    Standard_EXPORT GeomAbs_Shape Continuity() const Standard_OVERRIDE;
    //! 是否满足 C^N 连续
    Standard_EXPORT Standard_Boolean IsCN(const Standard_Integer N) const Standard_OVERRIDE;

    //! 仅计算点
    Standard_EXPORT void D0(const Standard_Real U, gp_Pnt& P) const Standard_OVERRIDE;
    //! 计算点 + 一阶导数
    Standard_EXPORT void D1(const Standard_Real U, gp_Pnt& P, gp_Vec& V1) const Standard_OVERRIDE;
    //! 计算点 + 一二阶导数
    Standard_EXPORT void D2(const Standard_Real U, gp_Pnt& P, gp_Vec& V1, gp_Vec& V2) const Standard_OVERRIDE;
    //! 计算点 + 一二三阶导数
    Standard_EXPORT void D3(const Standard_Real U, gp_Pnt& P, gp_Vec& V1, gp_Vec& V2, gp_Vec& V3) const Standard_OVERRIDE;
    //! N阶导数
    Standard_EXPORT gp_Vec DN(const Standard_Real U, const Standard_Integer N) const Standard_OVERRIDE;

    //! 对几何做变换
    Standard_EXPORT void Transform(const gp_Trsf& T) Standard_OVERRIDE;
    //! 对参数 U 的变换关系(简单返回 U)
    Standard_EXPORT virtual Standard_Real TransformedParameter(const Standard_Real U, const gp_Trsf& T) const Standard_OVERRIDE;
    //! 返回一个系数(简单返回 1.0)
    Standard_EXPORT virtual Standard_Real ParametricTransformation(const gp_Trsf& T) const Standard_OVERRIDE;

    //! 复制自身
    Standard_EXPORT Handle(Geom_Geometry) Copy() const Standard_OVERRIDE;

    //! 调试输出
    Standard_EXPORT  void DumpJson(Standard_OStream& theOStream, Standard_Integer theDepth = -1) const;

    //! RTTI 宏 (在头文件里)
    DEFINE_STANDARD_RTTIEXT(Geom_CurveOnSurface, Geom_BoundedCurve)

protected:
    //! 默认构造设为 protected
    Geom_CurveOnSurface() {}

private:
    Handle(Geom2d_BoundedCurve) myPCurve;   //!< 2D 曲线 (u,v)
    Handle(Geom_Surface) mySurface;  //!< 3D 曲面
};

#endif // _Geom_CurveOnSurface_HeaderFile