#include "CuttingRailsImpl.hxx"
#include <iostream>
#include <string>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stack>
#include <queue>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_TEdge.hxx>

#include <TopExp_Explorer.hxx>
#include <TColStd_Array1OfListOfInteger.hxx>
#include <TColStd_PackedMapOfInteger.hxx>

#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>

#include <GC_MakeCircle.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_TrimmedCurve.hxx>

#include <Geom_Plane.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>

#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gp_Lin.hxx>
#include <gp_Ax1.hxx>
#include <gp_Cylinder.hxx>

#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <ShapeFix_Wire.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_IntCS.hxx>
#include <NCollection_UBTreeFiller.hxx>
#include <BRepBuilderAPI_BndBoxTreeSelector.hxx>
#include <BndLib_Add3dCurve.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <ElSLib.hxx>

#include <TopExp.hxx>
#include <Adaptor3d_HCurve.hxx>

#include <math_Matrix.hxx>
#include <math_BullardGenerator.hxx>
#include <BRepClass_FaceExplorer.hxx>
#include <BRepClass_FClassifier.hxx>
#include <ElCLib.hxx>


const char* CuttingRailsImpl::BRepTools_Write(const char* theFileStr, void* theShapePtr)
{
	if (theFileStr == 0 || theShapePtr == 0)
	{
		return "Error: name or shape is null";
	}
	try {
		if (BRepTools::Write(*(TopoDS_Shape*)theShapePtr, theFileStr))
			return theFileStr;
		else
			return "Error: write failed";
			
	}
	catch (Standard_Failure const& anException)
	{
		return anException.GetMessageString();
	}
}

TopoDS_Shape CuttingRailsImpl::GetFirstSubShape(const TopoDS_Shape& shape) {
	// 首先检查输入是否有效
	if (shape.IsNull()) {
		return TopoDS_Shape();  // 返回空shape
	}

	TopoDS_Iterator exp(shape);

	// 检查是否有子形体
	if (!exp.More()) {
		return TopoDS_Shape();  // 没有子形体，返回空shape
	}

	return exp.Value();
}

Standard_Real CuttingRailsImpl::DistancePointToEdge(const gp_Pnt& point, const TopoDS_Edge& edge, const Standard_Integer nbPnts)
{
	Standard_Real f, l;
	Handle(Geom_Curve) gcurve = BRep_Tool::Curve(edge, f, l);
	// 离散采样
	double step = (l - f) / (nbPnts - 1);
	double param = f;
	gp_Pnt ptOnEdge;
	double minDist = 100000.0;
	for (size_t i = 0; i < nbPnts; i++)
	{
		param = f + step * i;
		gcurve->D0(param, ptOnEdge);
		Standard_Real dist = point.Distance(ptOnEdge);
		if (dist < minDist) {
			minDist = dist;
		}
	}

	return minDist;
}

// 计算点到Face边界的最小距离
Standard_Real CuttingRailsImpl::DistancePointToFaceBoundary(const gp_Pnt& point, const TopoDS_Face& face, double thinDist)
{
	Standard_Real minDist = std::numeric_limits<Standard_Real>::max();
	for (TopExp_Explorer expEdge(face, TopAbs_EDGE); expEdge.More(); expEdge.Next()) {
		TopoDS_Edge edge = TopoDS::Edge(expEdge.Current());

		Standard_Real f, l;
		Handle(Geom_Curve) gcurve = BRep_Tool::Curve(edge, f, l);
		if (gcurve.IsNull()) {
			return std::numeric_limits<Standard_Real>::max();
		}

		BRepAdaptor_Curve curve(edge);
		Standard_Real totalLen = GCPnts_AbscissaPoint::Length(curve); // edge曲线的长度
		int nbSamplesPerEdge = int(totalLen / thinDist) + 1;
		nbSamplesPerEdge = std::min(100, nbSamplesPerEdge);

		Standard_Real dist = DistancePointToEdge(point, edge, nbSamplesPerEdge);
		if (dist < minDist) {
			minDist = dist;
		}
	}
	return minDist;
}

double CuttingRailsImpl::FindMaxDistancePointOnFace(const TopoDS_Face& face, double thinDist, double deflection)
{
	// 1) 对Face进行网格剖分
	BRepMesh_IncrementalMesh mesher(face, deflection);

	// 2) 从面中提取网格(三角化)
	TopLoc_Location loc;
	Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, loc);
	if (triangulation.IsNull()) {
		std::cerr << "Error: Face Triangulation is null!\n";
		return -1.0;
	}

	// 3) 枚举网格顶点，计算到边界距离
	Standard_Real maxDist = 0.0;
	gp_Pnt bestPnt;
	int nbNode = triangulation->NbNodes();
	for (Standard_Integer i = 1; i <= triangulation->NbNodes(); ++i) {
		// 网格顶点在局部坐标 loc 下，需要变换回全局坐标
		gp_Pnt pnt = triangulation->Node(i).Transformed(loc.Transformation());

		// 计算此点到Face边界的最小距离
		Standard_Real dist = DistancePointToFaceBoundary(pnt, face, thinDist);
		if (dist > maxDist) {
			maxDist = dist;
			bestPnt = pnt;
		}
	}

	// 4) 返回找到的那个近似点
	return maxDist;
}


Standard_Boolean CuttingRailsImpl::FaceNormal(const TopoDS_Face& aF,
	const Standard_Real U,
	const Standard_Real V,
	gp_Dir& aDN)
{
	gp_Pnt aPnt;
	gp_Vec aD1U, aD1V, aN;
	Handle(Geom_Surface) aS;

	aS = BRep_Tool::Surface(aF);
	aS->D1(U, V, aPnt, aD1U, aD1V);
	aN = aD1U.Crossed(aD1V);
	if (aN.Magnitude() <= gp::Resolution())
		return Standard_False;

	aN.Normalize();
	aDN.SetXYZ(aN.XYZ());
	if (aF.Orientation() == TopAbs_REVERSED) {
		aDN.Reverse();
	}
	return Standard_True;
}

CuttingRailsImpl::ContinueEdge::ContinueEdge(const std::vector<TopoDS_Edge>& es, double length) : len(length)
{
	BRep_Builder BB;
	BB.MakeCompound(compedges);
	for (const auto& e : es) {

		BB.Add(compedges, e);
		TopoDS_Edge edge = TopoDS::Edge(e);
		if (BRep_Tool::Degenerated(edge)) {
			edges.push_back(edge);
		}
		else {
			edges.push_back(edge);
		}
	}
}

CuttingRailsImpl::ContinueEdge::ContinueEdge(const std::vector<TopoDS_Edge>& es, const TopoDS_Shape& ff, bool bOutLoop, double length, int index) : face(ff), faceIndex(index) {
	bOuterLoop = bOutLoop;
	degenerated = false;
	edgetType = EdgeType::NOTKNOWN;
	edges.insert(edges.begin(), es.begin(), es.end());
	len = length;

	BRep_Builder BB;
	BB.MakeCompound(compedges);
	for (const auto& e : edges) {
		BB.Add(compedges, e);
	}
}

CuttingRailsImpl::ContinueEdge::ContinueEdge(const std::vector<TopoDS_Edge>& es, const TopoDS_Shape& f, bool bOutLoop, int faceIndex) : face(f), faceIndex(faceIndex) {
	bOuterLoop = bOutLoop;
	degenerated = false;
	edgetType = EdgeType::NOTKNOWN;
	BRep_Builder BB;
	BB.MakeCompound(compedges);
	for (const auto& e : es) {

		BB.Add(compedges, e);
		TopoDS_Edge edge = TopoDS::Edge(e);
		if (BRep_Tool::Degenerated(edge)) {
			edges.push_back(edge);
		}
		else {
			edges.push_back(edge);
		}
	}
	if (edges.size() == 1 && BRep_Tool::Degenerated(edges[0])) {
		const TopoDS_Vertex& v1 = TopExp::FirstVertex(edges[0]);
		start = BRep_Tool::Pnt(v1);
		end = BRep_Tool::Pnt(v1);
		mid = BRep_Tool::Pnt(v1);
		len = 0.0;
		degenerated = true;
		bClosed = true;
	}
	else {
		const TopoDS_Vertex& v1 = TopExp::FirstVertex(edges[0], true);
		start.SetXYZ(BRep_Tool::Pnt(v1).XYZ());

		const TopoDS_Vertex& v2 = TopExp::LastVertex(edges[edges.size() - 1], true);
		end.SetXYZ(BRep_Tool::Pnt(v2).XYZ());

		bClosed = v1.IsSame(v2);

		double f, l;
		len = GetLength(increaseLen);
		double halfLen = len * 0.5;
		int midIndex = 0;
		for (int total = edges.size(); midIndex < total; midIndex++) {
			if (increaseLen[midIndex] >= halfLen) {
				break;
			}
		}
		const auto& edge = edges[midIndex];
		Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, f, l);
		//param = edge.Orientation() == face.Orientation() ? l : f;
		double param = (f + l) / 2;
		curve->D1(param, mid, midVec);
		if (edge.Orientation() == TopAbs_REVERSED) {
			midVec.Reverse();
		}
		midVec.Normalize();

		GetNormalOfMiddleEdge(midIndex);
	}
}

std::shared_ptr<CuttingRailsImpl::ContinueEdge> CuttingRailsImpl::ContinueEdge::Reversed() {

	std::vector<TopoDS_Edge> reversed_edges;
	reversed_edges.reserve(this->edges.size());

	// 反向遍历原vector，同时创建反向的edge
	for (auto it = this->edges.rbegin(); it != this->edges.rend(); ++it) {
		TopoDS_Edge reversed_edge = TopoDS::Edge(it->Reversed());
		reversed_edges.push_back(reversed_edge);
	}

	const auto& reversedContinueEdge = std::make_shared<ContinueEdge>(reversed_edges, this->face, this->bOuterLoop, this->faceIndex);
	reversedContinueEdge->mid = mid;
	reversedContinueEdge->midVec = midVec.Reversed();
	reversedContinueEdge->midNormal = midNormal.Reversed();
	return reversedContinueEdge;
}

TopoDS_Shape CuttingRailsImpl::ContinueEdge::GetCompound() {
	BRep_Builder BB;
	TopoDS_Compound compedges;
	BB.MakeCompound(compedges);
	for (const auto& e : edges) {
		BB.Add(compedges, e);
	}
	return compedges;

}

double CuttingRailsImpl::ContinueEdge::GetLength(std::vector<double>& increase) {
	double len = 0.0;
	for (size_t i = 0; i < edges.size(); i++)
	{
		GProp_GProps props;
		BRepGProp::LinearProperties(edges[i], props);
		len += props.Mass();
		increase.push_back(len);
	}
	return len;
}

void CuttingRailsImpl::ContinueEdge::GetNormalOfMiddleEdge(const int midIndex) {
	TopoDS_Edge edge;
	TopoDS_Face face = TopoDS::Face(this->face);
	Standard_Real f, l, param;
	gp_Pnt2d pt;
	Handle(Geom2d_Curve) curve2d;

	edge = this->edges[midIndex];
	curve2d = BRep_Tool::CurveOnSurface(edge, face, f, l);
	param = (f + l) / 2;
	curve2d->D0(param, pt);

	gp_Pnt aPnt;
	gp_Vec aD1U, aD1V, aN;
	Handle(Geom_Surface) aS;

	aS = BRep_Tool::Surface(face);
	aS->D1(pt.X(), pt.Y(), aPnt, aD1U, aD1V);
	aN = aD1U.Crossed(aD1V);

	aN.Normalize();
	midNormal.SetXYZ(aN.XYZ());
	if (face.Orientation() == TopAbs_REVERSED) {
		midNormal.Reverse();
	}
	return;
}

gp_Pnt2d CuttingRailsImpl::ContinueEdge::uvAtStart() {
	auto edge = edges[0];
	Standard_Real first, last;
	Handle(Geom2d_Curve) curve = BRep_Tool::CurveOnSurface(edge, TopoDS::Face(face), first, last);
	double param = edge.Orientation() == TopAbs_FORWARD ? first : last;

	gp_Pnt2d point;

	// 计算起点处的一阶导数
	curve->D0(param, point);

	return point;
}

bool CuttingRailsImpl::ContinueEdge::NormalAt(const gp_Pnt2d& pt2d, gp_Vec& normal) {
	gp_Pnt aPnt;
	gp_Vec aD1U, aD1V, aN, aDN;
	Handle(Geom_Surface) aS;

	aS = BRep_Tool::Surface(TopoDS::Face(face));
	aS->D1(pt2d.X(), pt2d.Y(), aPnt, aD1U, aD1V);
	aN = aD1U.Crossed(aD1V);
	if (aN.Magnitude() <= gp::Resolution()) {
		return false;
	}
	aN.Normalize();
	aDN.SetXYZ(aN.XYZ());
	if (face.Orientation() == TopAbs_REVERSED) {
		aDN.Reverse();
	}
	normal = aDN;
	return true;
}

gp_Pnt2d CuttingRailsImpl::ContinueEdge::uvAtEnd() {
	auto edge = edges[edges.size() - 1];
	Standard_Real first, last;
	Handle(Geom2d_Curve) curve = BRep_Tool::CurveOnSurface(edge, TopoDS::Face(face), first, last);

	gp_Pnt2d point;
	double param = edge.Orientation() == TopAbs_FORWARD ? last : first;

	// 计算终点处的一阶导数
	curve->D0(param, point);

	return point;
}

gp_Vec CuttingRailsImpl::ContinueEdge::dirAtStart() {
	auto edge = edges[0];
	Standard_Real first, last;
	Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);

	gp_Vec derivative;
	gp_Pnt point;

	double param = edge.Orientation() == TopAbs_FORWARD ? first : last;
	// 计算起点处的一阶导数
	curve->D1(param, point, derivative);

	if (derivative.Magnitude() <= gp::Resolution()) {
		double delta = (last - first) * 0.01;
		gp_Pnt p2;
		if (edge.Orientation() == TopAbs_FORWARD) {
			curve->D0(first + delta, p2);
			derivative = gp_Vec(point, p2);
		}
		else {
			curve->D0(last - delta, p2);
			derivative = gp_Vec(p2, point);
		}
	}
	// 如果边的方向是反向的，需要反转导数向量
	if (edge.Orientation() == TopAbs_REVERSED) {
		derivative.Reverse();
	}
	derivative.Normalize();
	return derivative;
}

gp_Vec CuttingRailsImpl::ContinueEdge::dirAtEnd() {
	auto edge = edges[edges.size() - 1];
	Standard_Real first, last;
	Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);

	gp_Vec derivative;
	gp_Pnt point;
	double param = edge.Orientation() == TopAbs_FORWARD ? last : first;

	// 计算起点处的一阶导数
	curve->D1(param, point, derivative);

	if (derivative.Magnitude() <= gp::Resolution()) {
		double delta = (last - first) * 0.01;
		gp_Pnt p2;
		if (edge.Orientation() == TopAbs_FORWARD) {
			curve->D0(last - delta, p2);
			derivative = gp_Vec(p2, point);
		}
		else {
			curve->D0(first + delta, p2);
			derivative = gp_Vec(point, p2);
		}
	}
	// 如果边的方向是反向的，需要反转导数向量
	if (edge.Orientation() == TopAbs_REVERSED) {
		derivative.Reverse();
	}
	derivative.Normalize();
	return derivative;
}

TopoDS_Vertex CuttingRailsImpl::ContinueEdge::GetStartVertex() {
	if (edges.empty()) {
		return TopoDS_Vertex();
	}
	return TopExp::FirstVertex(edges.front(), true);
}

TopoDS_Vertex CuttingRailsImpl::ContinueEdge::GetEndVertex() {
	if (edges.empty()) {
		return TopoDS_Vertex();
	}
	return TopExp::LastVertex(edges.back(), true);
}

bool CuttingRailsImpl::ContinueEdge::UpdateMidPt(const gp_Vec& midTangent, const gp_Pnt2d& midUV)
{
	Handle(Geom_Surface) surface = BRep_Tool::Surface(TopoDS::Face(face));
	if (surface.IsNull()) {
		return false;
	}
	gp_Vec aD1U, aD1V, aN;
	surface->D1(midUV.X(), midUV.Y(), mid, aD1U, aD1V);
	aN = aD1U.Crossed(aD1V);

	if (aN.SquareMagnitude() <= gp::Resolution()) {
		return false;
	}

	aN.Normalize();
	midNormal.SetXYZ(aN.XYZ());
	if (face.Orientation() == TopAbs_REVERSED) {
		midNormal.Reverse();
	}
	midVec = midTangent;
	return true;
}

TopoDS_Shape CuttingRailsImpl::showTopos(const TopTools_IndexedMapOfShape& shapes) {

	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);
	for (int i = 1; i <= shapes.Extent(); i++) {
		BB.Add(comp, shapes.FindKey(i));
	}
	return comp;
}

double CuttingRailsImpl::dist2VecVec(const point_t& p1, const point_t& pt) {
	double dx = p1.coords[0] - pt.coords[0];
	double dy = p1.coords[1] - pt.coords[1];
	double dz = p1.coords[2] - pt.coords[2];

	return dx * dx + dy * dy + dz * dz;
}

std::shared_ptr<KDTree> CuttingRailsImpl::MakeKDTree(const TopoDS_Shape& shape) {
	int i = 0;
	std::vector<point_t> vs;
	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next(), i++) {
		const TopoDS_Face& face = TopoDS::Face(exp.Current());
		TopTools_IndexedMapOfShape vMap;
		TopExp::MapShapes(face, TopAbs_VERTEX, vMap);
		for (Standard_Integer j = 1; j <= vMap.Extent(); j++) {
			const TopoDS_Vertex& v = TopoDS::Vertex(vMap(j));
			if (!v.IsNull()) {
				vs.push_back(point_t(v, i));
			}
		}
	}
	return std::make_shared<KDTree>(vs);
}

bool CuttingRailsImpl::IsConnectEdgeSmooth(const TopoDS_Shape& shape1, const TopoDS_Shape& shape2, const TopoDS_Vertex& v) {
	const auto& e1 = TopoDS::Edge(shape1);
	const auto& e2 = TopoDS::Edge(shape2);

	bool start1 = TopExp::FirstVertex(e1).IsSame(v);
	bool start2 = TopExp::FirstVertex(e2).IsSame(v);

	Standard_Real first, last;
	Handle(Geom_Curve) curve1 = BRep_Tool::Curve(e1, first, last);

	gp_Vec derivative1;
	gp_Pnt point;
	double param = start1 ? first : last;

	// 计算起点处的一阶导数
	curve1->D1(param, point, derivative1);

	// 如果边的方向是反向的，需要反转导数向量
	if (e1.Orientation() == TopAbs_REVERSED) {
		derivative1.Reverse();
	}
	derivative1.Normalize();

	Handle(Geom_Curve) curve2 = BRep_Tool::Curve(e2, first, last);

	gp_Vec derivative2;
	double param2 = start2 ? first : last;

	// 计算起点处的一阶导数
	curve2->D1(param2, point, derivative2);

	// 如果边的方向是反向的，需要反转导数向量
	if (e2.Orientation() == TopAbs_REVERSED) {
		derivative2.Reverse();
	}
	derivative2.Normalize();

	double dot = derivative1.Dot(derivative2);
	return dot > sqrt(3) / 2;
}

bool CuttingRailsImpl::IsPtBNearestPtA(const std::shared_ptr<KDTree>& tree, const point_t& ptA, const point_t& ptB, double tol) {
	if (ptB.attribute != ptA.attribute) {//查找P点不属于当前面的最近点.记为PP
		if (ptB.coords.empty() || ptA.coords.empty()) {
			// 错误处理
			return false;
		}
		double dist2 = dist2VecVec(ptB, ptA);
		if (dist2 < tol * tol) {//另一点和A小于给定容差
			const auto& ret2 = tree->nearest_points(ptB, 10);
			for (const auto& ele : ret2)
			{
				if (ele.attribute == ptA.attribute) {//另一点在当前面的最近点，
					if (abs(ele.coords[0] - ptA.coords[0]) < 1.0e-12 &&
						abs(ele.coords[1] - ptA.coords[1]) < 1.0e-12 &&
						abs(ele.coords[2] - ptA.coords[2]) < 1.0e-12) {
						return true;
					}
					else {
						return false;
					}
				}
			}

		}
		else {//PP和P大于给定容差
			return false;
		}
	}
	return false;
}

bool CuttingRailsImpl::GetMutualNearestPt(const std::shared_ptr<KDTree>& tree, const point_t& ptA, const pointVec& ptsVec, const double tol, point_t& ptNearest) {
	for (size_t k = 0; k < ptsVec.size(); k++)
	{
		if (ptA.attribute == ptsVec[k].attribute) {
			continue;
		}

		if (IsPtBNearestPtA(tree, ptA, ptsVec[k], tol)) {
			ptNearest = ptsVec[k];
			return true;
		}
		else {
			return false;
		}
	}
	return false;
}

gp_Vec CuttingRailsImpl::DirAtStart(const TopoDS_Edge& edge) {
	Standard_Real first, last;
	Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);


	gp_Pnt point;
	double delta = (last - first) * 0.03;
	double param1;
	double param2;
	if (edge.Orientation() == TopAbs_FORWARD) {
		param1 = first;
		param2 = first + delta;
	}
	else {
		param1 = last - delta;
		param2 = last;
	}
	// 计算起点处的一阶导数
	gp_Pnt p1, p2;
	curve->D0(param1, p1);
	curve->D0(param2, p2);
	gp_Vec derivative(p1, p2);
	// 如果边的方向是反向的，需要反转导数向量
	if (edge.Orientation() == TopAbs_REVERSED) {
		derivative.Reverse();
	}
	derivative.Normalize();
	return derivative;
}

bool CuttingRailsImpl::IsSmoothAtVertex(const TopTools_IndexedDataMapOfShapeListOfShape& VEMap, const TopoDS_Vertex& vertex) {
	const auto& edges = VEMap.FindFromKey(vertex);
	std::vector<gp_Dir> es;
	bool hasDegerated = false;
	double smoothDot = -cos(M_PI / 36);
	for (auto& ele : edges) {
		TopoDS_Edge edge = TopoDS::Edge(ele);
		if (!BRep_Tool::Degenerated(TopoDS::Edge(edge))) {
			if (TopExp::FirstVertex(edge, true).IsSame(vertex)) {
				es.push_back(DirAtStart(edge));
			}
			if (TopExp::LastVertex(edge, true).IsSame(vertex)) {
				es.push_back(DirAtStart(TopoDS::Edge(edge.Reversed())));
			}
		}
		else {
			hasDegerated = true;
		}
	}

	bool smooth = false;
	if (!hasDegerated && es.size() == 2) {
		smooth = es[0].Dot(es[1]) < smoothDot; // 30度
	}
	return smooth;
}

bool CuttingRailsImpl::MakeContinueEdges(const std::shared_ptr<KDTree>& tree, const TopoDS_Shape& shape,
	std::vector<ContinueEdgePtr>& wireContinueEdges, TopTools_IndexedMapOfShape& softVs, double tol) {
	wireContinueEdges.clear();
	int i = 0;
	double smoothDot = -cos(M_PI / 36);//- sqrt(3) / 2;
	TopTools_IndexedDataMapOfShapeListOfShape VEMap;
	TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, VEMap);
	TopTools_IndexedMapOfShape allRealVertexMap;
	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next(), i++) {
		TopoDS_Shape face = exp.Current();
		//begin collect real vretices
		TopTools_IndexedMapOfShape vMap;
		TopTools_IndexedMapOfShape realVertexs;

		TopExp::MapShapes(face, TopAbs_VERTEX, vMap);
		for (Standard_Integer j = 1; j <= vMap.Extent(); j++)
		{
			const TopoDS_Vertex& v = TopoDS::Vertex(vMap(j));
			if (v.IsNull()) {
				continue;
			}

			if (softVs.Contains(v)) {
				continue;
			}

			if (allRealVertexMap.Contains(v)) {
				realVertexs.Add(v);
				continue;
			}
			auto& vEdges = VEMap.FindFromKey(v);
			bool hasDegenerateEdge = false;
			for (auto& e : vEdges) {
				if (BRep_Tool::Degenerated(TopoDS::Edge(e))) {
					realVertexs.Add(v);
					hasDegenerateEdge = true;
					break;
				}
			}
			if (hasDegenerateEdge) {
				continue;
			}

			//判断点是不是real vertex
			//双向匹配: faceA上的一点和faceB上的一点互为最近点
			point_t ptInfoA(v, i);
			TopTools_IndexedMapOfShape matchedVertices;
			matchedVertices.Add(v);
			std::set<int> anotherFaceIndexSet;// 确保A点在另一个面上双向匹配  只匹配一个点
			const auto& ret = tree->nearest_points(ptInfoA, 10);
			for (size_t k = 0; k < ret.size(); k++)
			{
				if (anotherFaceIndexSet.find(ret[k].attribute) == anotherFaceIndexSet.end()) {
					if (IsPtBNearestPtA(tree, ptInfoA, ret[k], tol)) {
						anotherFaceIndexSet.insert(ret[k].attribute);
						matchedVertices.Add(ret[k].vertex);
					}
				}
			}

			bool isRealVertex = matchedVertices.Extent() >= 2 ||
				(matchedVertices.Extent() == 1 && vEdges.Size() > 2);

			if (isRealVertex) {
				if (matchedVertices.Contains(v)) {
					realVertexs.Add(v);
				}
				for (int i = 1; i <= matchedVertices.Extent(); i++) {
					allRealVertexMap.Add(matchedVertices(i));
				}
			}
		}

		//auto realVsShape = showTopos(realVertexs);
		//end collect real vretices

		//begin splitwires by realvertex
		TopoDS_Wire outerWire = BRepTools::OuterWire(TopoDS::Face(face));
		for (TopExp_Explorer exp(face, TopAbs_WIRE); exp.More(); exp.Next()) {
			const TopoDS_Wire& wire = TopoDS::Wire(exp.Current());
			std::vector<ContinueEdgePtr> continueEdges;
			TopExp_Explorer wireExplorer(wire, TopAbs_EDGE);
			int index = -1;
			//make edges form wire
			std::vector<TopoDS_Edge> edges;
			for (; wireExplorer.More(); wireExplorer.Next()) {
				const TopoDS_Edge& edge = TopoDS::Edge(wireExplorer.Current());

				if (index == -1) {
					TopoDS_Vertex v1 = TopExp::FirstVertex(edge, true);
					bool bV1RealVertex = realVertexs.Contains(v1);
					if (bV1RealVertex) {
						index = edges.size();
					}
				}
				edges.push_back(edge);
			}

			int size = edges.size();
			int startIndex = index == -1 ? 0 : index;
			std::vector<TopoDS_Edge> continues;
			std::vector<int> eIndex;
			for (int ii = startIndex, len = startIndex + size; ii < len; ii++) {
				int id = ii >= size ? ii - size : ii;
				const TopoDS_Edge& edge = edges[ii >= size ? ii - size : ii];
				bool bDegenerated = BRep_Tool::Degenerated(edge);
				if (bDegenerated) {
					if (!continues.empty()) {
						wireContinueEdges.push_back(std::make_shared<ContinueEdge>(continues, face, wire.IsEqual(outerWire), i));
						continues.clear();
						eIndex.clear();
					}
					eIndex.push_back(id);
					wireContinueEdges.push_back(std::make_shared<ContinueEdge>(std::vector<TopoDS_Edge>{edge}, face, wire.IsEqual(outerWire), i));
					eIndex.clear();
				}
				else {
					// 检查端点的vertexCount
					continues.push_back(edge);
					eIndex.push_back(id);
					bool bV2RealVertex = realVertexs.Contains(TopExp::LastVertex(edge, true));
					if (bV2RealVertex) {
						wireContinueEdges.push_back(std::make_shared<ContinueEdge>(continues, face, wire.IsEqual(outerWire), i));
						continues.clear();
						eIndex.clear();
					}
				}
			}

			if (index == -1) {
				// loop 不存在非soft点，按理不应该发生
				wireContinueEdges.push_back(std::make_shared<ContinueEdge>(continues, face, wire.IsEqual(outerWire), i));
			}
			else if (!continues.empty()) {
				wireContinueEdges.push_back(std::make_shared<ContinueEdge>(continues, face, wire.IsEqual(outerWire), i));
			}

			//end splitwires by realvertex
		}
	}
	return true;
}

void CuttingRailsImpl::GetEdgeCurveTangent(const TopoDS_Shape& edge, const Handle(Geom_Curve)& curve, double param, gp_Vec& tangent)
{
	gp_Pnt pt;
	curve->D1(param, pt, tangent);
	if (edge.Orientation() == TopAbs_REVERSED) {
		tangent.Reverse();
	}
	tangent.Normalize();
}

size_t CuttingRailsImpl::GetVertexCount(const std::vector<std::pair<gp_Pnt, int>>& vertexCount, const gp_Pnt& pt, double tol) {
	for (size_t i = 0; i < vertexCount.size(); i++) {
		if (pt.Distance(vertexCount[i].first) < tol) {
			return vertexCount[i].second;
		}
	}
	return 0;
}

size_t CuttingRailsImpl::GetVertexCount(const std::vector<std::pair<gp_Pnt, std::vector<TopoDS_Edge>>>& vertexCount, const gp_Pnt& pt) {
	for (size_t i = 0; i < vertexCount.size(); i++) {
		if (pt.Distance(vertexCount[i].first) < Precision::Confusion()) {
			return vertexCount[i].second.size();
		}
	}
	return 0;
}

void CuttingRailsImpl::PrintFaceEdge(const std::vector<ContinueEdgePtr>& wireContinueEdges, int fIndex) {
	std::string fileName = "E:/temp_test/";
	int eIndex = 1;
	for (auto& continues : wireContinueEdges) {
		auto& edges = continues->edges;
		for (auto& edge : edges) {
			std::string path = fileName + "continue_f" + std::to_string(fIndex) + "e" + std::to_string(eIndex) + ".brep";
			BRepTools::Write(edge, path.c_str());
			eIndex++;
		}
	}
}

double CuttingRailsImpl::PointDistanceToCompCurve(const gp_Pnt& p, const ContinueEdgePtr& compCurve, gp_Vec& refMidVec, gp_Pnt2d& refMidUV)
{
	// compCurve 可能包含若干子曲线(对应原先的edge)
	// 使用 compCurve.Curve(i) 取第 i 条子曲线
	double minDist = 100000;
	const size_t nbc = compCurve->edges.size();
	const TopoDS_Face compFace = TopoDS::Face(compCurve->face);
	gp_Pnt mid;

	for (int i = 0; i < nbc; ++i) {

		// 取得它在 BRep 中对应的 Edge 和有效参数范围 [f, l]
		const TopoDS_Edge& subEdge = TopoDS::Edge(compCurve->edges[i]);
		Standard_Real fParam;
		Standard_Real lParam;
		Standard_Real fpParam;
		Standard_Real lpParam;

		// 从 Edge 拿到底层 Geom_Curve(可能是 Line/Circle/BSpline等)
		// 注意 BRep_Tool::Curve 返回 Handle(Geom_Curve)
		Handle(Geom_Curve) hCurve = BRep_Tool::Curve(subEdge, fParam, lParam);

		Handle(Geom2d_Curve) pCurve = BRep_Tool::CurveOnSurface(subEdge, compFace, fpParam, lpParam);
		if (hCurve.IsNull() || pCurve.IsNull()) {
			continue;
		}
		// 用 GeomAPI_ProjectPointOnCurve 把点投影到这条曲线
		// 并获取最小距离
		GeomAPI_ProjectPointOnCurve projector;
		projector.Init(p, hCurve, fParam, lParam);
		if (projector.NbPoints() > 0) {
			double projectMinDist = projector.LowerDistance();
			if (minDist > projectMinDist) {
				minDist = projectMinDist;
				double minParam = projector.LowerDistanceParameter();
				pCurve->D0(minParam, refMidUV);
				GetEdgeCurveTangent(subEdge, hCurve, minParam, refMidVec);
			}
		}

		TopoDS_Vertex v1 = TopExp::FirstVertex(subEdge);
		TopoDS_Vertex v2 = TopExp::LastVertex(subEdge);

		double d1 = BRep_Tool::Pnt(v1).Distance(p);
		double d2 = BRep_Tool::Pnt(v2).Distance(p);
		if (d1 < minDist) {
			minDist = d1;
			pCurve->D0(fpParam, refMidUV);
			GetEdgeCurveTangent(subEdge, hCurve, fParam, refMidVec);
		}

		if (d2 < minDist) {
			minDist = d2;
			pCurve->D0(lpParam, refMidUV);
			GetEdgeCurveTangent(subEdge, hCurve, lParam, refMidVec);
		}
	}

	// 如果所有子曲线都查完，没有距离<=tol的情况，返回 false
	return minDist;
}

int CuttingRailsImpl::NbChildren(const TopoDS_Shape& shape)
{
	TopoDS_Iterator it(shape);
	int count = 0;
	for (; it.More(); it.Next()) {
		++count;
	}
	return count;
}

TopoDS_Compound CuttingRailsImpl::GetUnMatchedEdges(const std::vector<ContinueEdgePtr>& CCUnMatched, const double tol) {
	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);

	for (const auto coedges : CCUnMatched)
	{
		if (coedges->len < 0.01 * tol) {
			continue;//小于容差的边不进行匹配
		}
		BB.Add(comp, coedges->GetCompound());
	}

	return comp;
}

int CuttingRailsImpl::MatchBFromA(const std::shared_ptr<KDTree>& tree,
	const VertexToContinueMap& vertexToContinueMap,
	const size_t ce1Index,
	const std::vector<ContinueEdgePtr>& continueEdges,
	double tol,
	double midTol, bool& sameDir, gp_Vec& matchedMidVec, gp_Pnt2d& matchedMidUV) {
	auto& ce1 = continueEdges[ce1Index];
	auto& v1 = ce1->GetStartVertex();
	auto& v2 = ce1->GetEndVertex();

	int ce1FaceIndex = ce1->faceIndex;

	point_t ptInfo1(v1, ce1FaceIndex);
	point_t ptInfo2(v2, ce1FaceIndex);

	const std::vector<point_t>& ret1 = tree->nearest_points(ptInfo1, 10);
	const std::vector<point_t>& ret2 = tree->nearest_points(ptInfo2, 10);

	std::unordered_set<size_t> resultSet1, resultSet2;

	for (const auto& pt : ret1) {
		if (vertexToContinueMap.Contains(pt.vertex)) {
			const std::vector<size_t>& ids = vertexToContinueMap.FindFromKey(pt.vertex);
			resultSet1.insert(ids.begin(), ids.end());  // 插入并自动去重
		}
	}

	for (const auto& pt : ret2) {
		if (vertexToContinueMap.Contains(pt.vertex)) {
			const std::vector<size_t>& ids = vertexToContinueMap.FindFromKey(pt.vertex);
			resultSet2.insert(ids.begin(), ids.end());  // 插入并自动去重
		}
	}


	std::vector<size_t> result;
	for (const auto& p : resultSet1) {
		auto it = std::find_if(resultSet2.begin(), resultSet2.end(), [&](const size_t& q) {
			return p == q;  // 比较的是内容
			});
		if (it != resultSet2.end()) {
			result.push_back(p);
		}
	}

	// 肯定存在一条当前ce1
	if (result.size() >= 2) {
		int matchedIndex = -1;
		double matchMinDistinct = 1000000;
		gp_Vec ce2MidVec;
		gp_Pnt2d ce2MidUV;
		for (auto& j : result) {
			size_t ce2Index = j;
			auto& ce2 = continueEdges[ce2Index];
			if (ce1Index == ce2Index) {
				continue;
			}
			if (ce1->faceIndex == ce2->faceIndex) {
				if (ce1->edges.size() == 1 && ce2->edges.size() == 1) {
					// 缝合边判断
					if (ce1->edges[0].IsSame(ce2->edges[0])) {
						return ce2Index;
					}
				}
				continue;
			}

			double endPntMatch1 = (ce1->start.SquareDistance(ce2->start) + ce1->end.SquareDistance(ce2->end));
			double endPntMatch2 = (ce1->start.SquareDistance(ce2->end) + ce1->end.SquareDistance(ce2->start));
			bool endPntMatch = endPntMatch1 < 2 * tol * tol || endPntMatch2 < 2 * tol * tol;
			if (endPntMatch) {
				//找最近的匹配边
				double midDist1 = PointDistanceToCompCurve(ce1->mid, ce2, ce2MidVec, ce2MidUV);

				if (midDist1 < midTol) {
					// 直接计算新的匹配距离
					double addMid = 0.1 * midDist1;
					double newMatchDist = std::min(endPntMatch1 + addMid, endPntMatch2 + addMid);

					if (matchMinDistinct > newMatchDist) {
						matchMinDistinct = newMatchDist;

						sameDir = endPntMatch1 < 2 * tol * tol;
						if (endPntMatch1 < 2 * tol * tol && endPntMatch2 < 2 * tol * tol)
						{
							sameDir = ce1->dirAtStart().Dot(ce2->dirAtStart()) > 0.0 ? true : false;
						}
						matchedIndex = ce2Index;
						matchedMidVec = ce2MidVec;
						matchedMidUV = ce2MidUV;
					}
				}
			}
		}
		return matchedIndex;
	}
	else {
		return -1;
	}
}

bool CuttingRailsImpl::MatchCoedges(const std::shared_ptr<KDTree>& tree, const std::vector<ContinueEdgePtr>& continueEdges, CoedgeInfo& CCMatchedMap, std::vector<ContinueEdgePtr>& CCUnMatched, double tol, TopoDS_Compound& unMatchedEdges) {
	CCMatchedMap.clear();
	CCUnMatched.clear();

	// 连续边匹配
	std::set<int> usedFlag;
	double midTol = std::min(tol * 10, 0.5);
	VertexToContinueMap vertexToContinueMap;
	// 建立 vertex 到 ContinueEdge 的映射, map内不含退化的连续边和小于容差的连续边
	for (size_t i = 0, iLen = continueEdges.size(); i < iLen; i++) {
		auto& ce = continueEdges[i];
		if (ce->degenerated) {
			CCMatchedMap.push_back({ std::vector< ContinueEdgePtr>{ce, ce}, false });
			continue;
		}

		if (ce->len < 0.01 * tol) {
			CCUnMatched.push_back(ce);
			continue;//小与容差的边不进行匹配
		}

		auto& start = ce->GetStartVertex();
		if (vertexToContinueMap.Contains(start)) {
			vertexToContinueMap.ChangeFromKey(start).push_back(i);
		}
		else {
			vertexToContinueMap.Add(start, { i });
		}

		auto& end = ce->GetEndVertex();
		if (vertexToContinueMap.Contains(end)) {
			vertexToContinueMap.ChangeFromKey(end).push_back(i);
		}
		else {
			vertexToContinueMap.Add(end, { i });
		}
	}

	for (int i = 0, iLen = continueEdges.size(); i < iLen; i++) {
		auto& ce = continueEdges[i];
		if (usedFlag.find(i) != usedFlag.end()) {
			// 如果已经使用过，跳过
			continue;
		}

		if (ce->degenerated) {
			continue;
		}

		if (ce->len < 0.01 * tol) {
			continue;//小与容差的边不进行匹配
		}

		bool sameDir = false;
		gp_Pnt2d ce1MidUV, ce2MidUV;
		gp_Vec ce1MidVec, ce2MidVec;
		int matchedIndex = MatchBFromA(tree, vertexToContinueMap, i, continueEdges, tol, midTol, sameDir, ce2MidVec, ce2MidUV);
		if (matchedIndex == -1 || usedFlag.find(matchedIndex) != usedFlag.end()) {
			CCUnMatched.push_back(continueEdges[i]);
			usedFlag.insert(i);
		}
		else {
			int matchedIndex2 = MatchBFromA(tree, vertexToContinueMap, matchedIndex, continueEdges, tol, midTol, sameDir, ce1MidVec, ce1MidUV);
			if (matchedIndex2 == -1 || matchedIndex2 != i) {
				// 双向匹配失败
				CCUnMatched.push_back(continueEdges[i]);
			}
			else {
				CCMatchedMap.push_back({ { continueEdges[i], continueEdges[matchedIndex] }, sameDir });
				continueEdges[matchedIndex]->UpdateMidPt(ce2MidVec, ce2MidUV);
				usedFlag.insert(matchedIndex);
			}
			usedFlag.insert(i);
		}
	}

	unMatchedEdges = GetUnMatchedEdges(CCUnMatched, tol);
	// BRepTools_Write("E:/unmatched.brep", (void*)&unMatchedEdges);
	//std::cout << "Unmatched Edges Count: " << NbChildren(unMatchedEdges) << std::endl;

	return true;
}

bool CuttingRailsImpl::CanMerge(const ContinueEdgePtr& continue1, const ContinueEdgePtr& continue2, bool checkEnd) {
	if (checkEnd) {
		return TopExp::LastVertex(continue1->edges.back(), true).IsSame(TopExp::FirstVertex(continue2->edges.front(), true));
	}
	else {
		return TopExp::FirstVertex(continue1->edges.front(), true).IsSame(TopExp::LastVertex(continue2->edges.back(), true));
	}
}

const CuttingRailsImpl::ContinueEdgePtr CuttingRailsImpl::MergeEdges(const std::vector<ContinueEdgePtr>& continuesToMerge) {
	if (continuesToMerge.size() == 1) {
		return continuesToMerge[0];
	}

	std::vector<TopoDS_Edge> edges;
	double totalLen = 0;
	std::vector<double> increase;
	for (auto& singleContinue : continuesToMerge) {
		edges.insert(edges.end(), singleContinue->edges.begin(), singleContinue->edges.end());
		totalLen += singleContinue->len;
		increase.push_back(totalLen);
	}

	ContinueEdgePtr newContinue = std::make_shared<ContinueEdge>(edges, continuesToMerge[0]->face, continuesToMerge[0]->bOuterLoop, totalLen, continuesToMerge[0]->faceIndex);
	newContinue->start = continuesToMerge.front()->start;
	newContinue->end = continuesToMerge.back()->end;

	double f, l;
	int size = edges.size();
	int middleIndex = 0;
	double halfLen = totalLen * 0.5;
	for (int count = continuesToMerge.size(); middleIndex < count; middleIndex++) {
		if (increase[middleIndex] >= halfLen) {
			break;
		}
	}

	newContinue->mid = continuesToMerge[middleIndex]->mid;
	newContinue->midVec = continuesToMerge[middleIndex]->midVec;
	newContinue->midNormal = continuesToMerge[middleIndex]->midNormal;
	return newContinue;
}

std::vector<CuttingRailsImpl::ContinueEdgePtr> CuttingRailsImpl::MergeGroup(const std::vector<ContinueEdgePtr>& group) {
	if (group.size() <= 1)
		return group;

	std::vector<ContinueEdgePtr> result;
	std::vector<bool> merged(group.size(), false);

	for (size_t i = 0; i < group.size(); i++) {
		if (merged[i])
			continue;

		std::vector<ContinueEdgePtr> continuesToMerge;
		continuesToMerge.push_back(group[i]);
		std::queue<size_t> toMergeQueue;
		toMergeQueue.push(i);
		merged[i] = true;

		// 后向查找
		while (!toMergeQueue.empty()) {
			auto& curContinueIndex = toMergeQueue.front();
			toMergeQueue.pop();
			const auto& curContinue = group[curContinueIndex];
			for (size_t j = 0; j < group.size(); j++) {
				if (curContinueIndex == j || merged[j])
					continue;
				if (CanMerge(curContinue, group[j], true)) {
					continuesToMerge.push_back(group[j]);
					merged[j] = true;
					toMergeQueue.push(j);
					break;
				}
			}
		}

		// 前向查找
		toMergeQueue.push(i);
		while (!toMergeQueue.empty()) {
			auto& curContinueIndex = toMergeQueue.front();
			toMergeQueue.pop();
			const auto& curContinue = group[curContinueIndex];
			for (size_t j = 0; j < group.size(); j++) {
				if (curContinueIndex == j || merged[j])
					continue;
				if (CanMerge(curContinue, group[j], false)) {
					continuesToMerge.emplace(continuesToMerge.begin(), group[j]);
					merged[j] = true;
					toMergeQueue.push(j);
				}
			}
		}


		result.push_back(MergeEdges(continuesToMerge));
	}

	return result;
}

const std::vector<CuttingRailsImpl::ContinueEdgePtr> CuttingRailsImpl::MergeContinuousEdges(const std::vector<ContinueEdgePtr>& continues) {
	if (continues.empty() || continues.size() == 1)
		return continues;

	std::vector<ContinueEdgePtr> result;
	// 按face和wire分组
	FaceToContinuesMap groupedContinues;
	for (const auto& singleContinue : continues) {
		if (groupedContinues.Contains(singleContinue->face)) {
			groupedContinues.ChangeFromKey(singleContinue->face).push_back(singleContinue);
		}
		else {
			groupedContinues.Add(singleContinue->face, std::vector<ContinueEdgePtr>{ singleContinue });
		}
	}

	// 处理每个组
	for (const auto& group : groupedContinues) {
		auto mergedEdges = MergeGroup(group);
		result.insert(result.end(), mergedEdges.begin(), mergedEdges.end());
	}

	return result;
}

bool CuttingRailsImpl::MatchMergedEdges(const std::shared_ptr<KDTree>& tree, const std::vector<ContinueEdgePtr>& FirstUnMatched,
	CoedgeInfo& CCMatchedMap, std::vector<ContinueEdgePtr>& LastUnMatched, double tol, TopoDS_Compound& unMatchedEdges) {
	LastUnMatched.clear();
	// 此时, 第二次边匹配后的结果会全部更新到CCMatchedMap, 包括匹配成功、未匹配成功、退化、短边
	const std::vector<ContinueEdgePtr> merged = MergeContinuousEdges(FirstUnMatched);

	CoedgeInfo secondMatchedMap;
	MatchCoedges(tree, merged, secondMatchedMap, LastUnMatched, tol, unMatchedEdges);
	if (!secondMatchedMap.empty()) {
		CCMatchedMap.insert(CCMatchedMap.end(), secondMatchedMap.begin(), secondMatchedMap.end());
	}

	return true;
}

double CuttingRailsImpl::CDistance(const TopoDS_Face& face) {
	// 创建线圈探索器
	TopExp_Explorer wireExp;
	std::vector<TopoDS_Wire> wires;

	// 收集face中的所有线圈
	for (wireExp.Init(face, TopAbs_WIRE); wireExp.More(); wireExp.Next()) {
		wires.push_back(TopoDS::Wire(wireExp.Current()));
	}

	// 如果线圈数量不是2，返回-1表示错误
	if (wires.size() != 2) return -1;

	std::vector<gp_Pnt> points1, points2;
	// 对每个线圈进行离散化
	for (int wireIndex = 0; wireIndex < 2; wireIndex++) {
		const TopoDS_Wire& wire = wires[wireIndex];
		std::vector<gp_Pnt>& points = (wireIndex == 0) ? points1 : points2;
		// 遍历线圈中的边
		for (TopExp_Explorer edgeExp(wire, TopAbs_EDGE); edgeExp.More(); edgeExp.Next()) {
			TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());
			// 获取边的几何信息
			Standard_Real first, last;
			Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);
			// 计算离散点的数量（这里使用固定数量，也可以根据曲线长度动态确定）
			int nbPoints = 50;
			Standard_Real delta = (last - first) / (nbPoints - 1);
			// 在边上生成离散点
			for (int i = 0; i < nbPoints; i++) {
				Standard_Real param = first + i * delta;
				gp_Pnt point;
				curve->D0(param, point);
				points.push_back(point);
			}
		}
	}

	// 计算从points1到points2的倒角距离之和
	double sum1 = 0.0;
	for (const gp_Pnt& p1 : points1) {
		double minDist = std::numeric_limits<double>::max();
		for (const gp_Pnt& p2 : points2) {
			double dist = p1.Distance(p2);
			minDist = std::min(minDist, dist);
		}
		sum1 += minDist;
	}

	// 计算从points2到points1的倒角距离之和  
	double sum2 = 0.0;
	for (const gp_Pnt& p2 : points2) {
		double minDist = std::numeric_limits<double>::max();
		for (const gp_Pnt& p1 : points1) {
			double dist = p2.Distance(p1);
			minDist = std::min(minDist, dist);
		}
		sum2 += minDist;
	}

	// 倒角距离是两个和的平均值
	return std::min(sum1 / points1.size(), sum2 / points2.size());
}

double CuttingRailsImpl::GetMaxLength(const TopoDS_Face& face, const std::vector<std::pair<gp_Pnt, int>>& vertexCount) {
	double maxLength = 0.0;
	double currentLength = 0.0;

	TopExp_Explorer wireExp(face, TopAbs_WIRE);
	if (!wireExp.More()) {
		return -1;  // 没有wire
	}

	const TopoDS_Wire& wire = TopoDS::Wire(wireExp.Current());

	// 使用BRepTools_WireExplorer保证edge的连续性
	TopExp_Explorer wireExplorer(wire, TopAbs_EDGE);
	int index = -1;
	std::vector<TopoDS_Edge> edges;
	for (; wireExplorer.More(); wireExplorer.Next()) {
		const TopoDS_Edge& edge = TopoDS::Edge(wireExplorer.Current());

		if (index == -1) {
			TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
			if (GetVertexCount(vertexCount, BRep_Tool::Pnt(v1), BRep_Tool::Tolerance(v1)) > 2) {
				index = edges.size();
			}
		}
		edges.push_back(edge);
	}

	int size = edges.size();
	int startIndex = index == -1 ? 0 : index;
	for (int i = startIndex, len = startIndex + size; i < len; i++) {
		const TopoDS_Edge& edge = edges[i >= size ? i - size : i];
		TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
		TopoDS_Vertex v2 = TopExp::LastVertex(edge);

		// 检查端点的vertexCount
		bool v2Count2 = GetVertexCount(vertexCount, BRep_Tool::Pnt(v2), BRep_Tool::Tolerance(v1)) == 2;

		// 计算edge长度
		GProp_GProps props;
		BRepGProp::LinearProperties(edge, props);
		double length = props.Mass();

		if (!v2Count2) {
			currentLength += length;
			if (currentLength > maxLength) {
				maxLength = currentLength;
			}
			currentLength = 0;
		}
		else {
			currentLength += length;
		}
	}
	return index == -1 ? currentLength : maxLength;
}

double CuttingRailsImpl::GetMaxLength(const TopoDS_Face& face, const std::vector<std::pair<gp_Pnt, std::vector<TopoDS_Edge>>>& vertexCount) {
	double maxLength = 0.0;
	double currentLength = 0.0;

	TopExp_Explorer wireExp(face, TopAbs_WIRE);
	if (!wireExp.More()) {
		return -1;  // 没有wire
	}

	const TopoDS_Wire& wire = TopoDS::Wire(wireExp.Current());

	// 使用BRepTools_WireExplorer保证edge的连续性
	TopExp_Explorer wireExplorer(wire, TopAbs_EDGE);
	int index = -1;
	std::vector<TopoDS_Edge> edges;
	for (; wireExplorer.More(); wireExplorer.Next()) {
		const TopoDS_Edge& edge = TopoDS::Edge(wireExplorer.Current());

		if (index == -1) {
			TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
			if (GetVertexCount(vertexCount, BRep_Tool::Pnt(v1)) > 2) {
				index = edges.size();
			}
		}
		edges.push_back(edge);
	}

	int size = edges.size();
	int startIndex = index == -1 ? 0 : index;
	for (int i = startIndex, len = startIndex + size; i < len; i++) {
		const TopoDS_Edge& edge = edges[i >= size ? i - size : i];
		TopoDS_Vertex v1 = TopExp::FirstVertex(edge);
		TopoDS_Vertex v2 = TopExp::LastVertex(edge);

		// 检查端点的vertexCount
		bool v2Count2 = GetVertexCount(vertexCount, BRep_Tool::Pnt(v2)) <= 2;

		// 计算edge长度
		GProp_GProps props;
		BRepGProp::LinearProperties(edge, props);
		double length = props.Mass();

		if (!v2Count2) {
			currentLength += length;
			if (currentLength > maxLength) {
				maxLength = currentLength;
			}
			currentLength = 0;
		}
		else {
			currentLength += length;
		}
	}
	return index == -1 ? currentLength : maxLength;
}

double CuttingRailsImpl::GetPerimeter(const TopoDS_Shape& comp) {

	double tol = 0.0;

	for (TopExp_Explorer expFace(comp, TopAbs_FACE); expFace.More(); expFace.Next())
	{

		for (TopExp_Explorer exp(expFace.Current(), TopAbs_EDGE); exp.More(); exp.Next()) {

			if (BRep_Tool::Degenerated(TopoDS::Edge(exp.Current()))) {
				continue;
			}
			GProp_GProps props;
			BRepGProp::LinearProperties(exp.Current(), props);
			tol += props.Mass();
		}
	}

	return tol;
}

std::vector<CuttingRailsImpl::ContinueEdgePtr> CuttingRailsImpl::GetCandidateRails(const TopoDS_Shape& outBoundFaces, const Face2CoedgesMap& faceCoedgesMap, bool needCuttingFace) {
	TopTools_IndexedMapOfShape map;
	for (TopExp_Explorer exp(outBoundFaces, TopAbs_FACE); exp.More(); exp.Next()) {
		map.Add(exp.Current());
	}
	std::vector<ContinueEdgePtr> candidatesRails;

	for (int i = 1; i <= map.Extent(); i++) {
		auto& curFace = map.FindKey(i);
		std::vector<ContinueEdgePtr> candidates;
		if (faceCoedgesMap.Contains(curFace)) {
			auto& continues = *(faceCoedgesMap.FindFromKey(curFace));
			for (auto& pair : continues) {
				auto& relate = pair.second;
				if (pair.second != nullptr && map.Contains(pair.second->face)) {
					continue;
				}
				else {
					candidates.push_back(pair.first);
				}
			}
		}
		candidatesRails.insert(candidatesRails.end(), candidates.begin(), candidates.end());
	}
	return candidatesRails;
}

bool CuttingRailsImpl::Contains(const VertexToContinueIndexMap& vertexMap, const gp_Pnt& point, const double squareTol) {
	for (auto& vertex : vertexMap) {
		if (point.SquareDistance(vertex.first) <= squareTol) {
			return true;
		}
	}
	return false;
}

std::vector<int> CuttingRailsImpl::GetValue(const VertexToContinueIndexMap& vertexMap, const gp_Pnt& point, const double squareTol) {
	std::vector<std::pair<int, double>> candidates;
	for (auto& vertex : vertexMap) {
		double sq = point.SquareDistance(vertex.first);
		if (sq <= squareTol) {
			candidates.push_back({ vertex.second, sq });
		}
	}
	std::sort(candidates.begin(), candidates.end(),
		[](const auto& a, const auto& b) {
			return a.second < b.second;
		});
	std::vector<int> result;
	result.reserve(candidates.size());  // 提前分配内存
	for (const auto& p : candidates) {
		result.push_back(p.first);
	}
	return result;
}

CuttingRailsImpl::VertexToContinueIndexMap CuttingRailsImpl::BuildVertexToContinueIndexMap(const std::vector<ContinueEdgePtr>& continues) {
	VertexToContinueIndexMap vertexMap;
	for (int i = 0, iLen = continues.size(); i < iLen; i++) {
		TopoDS_Vertex v1 = continues[i]->GetStartVertex();
		gp_Pnt p1 = BRep_Tool::Pnt(v1);
		vertexMap.push_back({ p1, i });
	}

	return vertexMap;
}

// 查找下一条边
int CuttingRailsImpl::FindNextEdge(const gp_Pnt& point,
	const VertexToContinueIndexMap& vertexMap,
	const std::set<int>& usedEdges,
	double squareTol) {
	// 查找与当前点相连的所有边
	if (!Contains(vertexMap, point, squareTol)) {
		return -1;
	}

	// 在所有相连的边中找到未使用的边
	auto& nextCandidates = GetValue(vertexMap, point, squareTol);
	for (auto& candidate : nextCandidates) {
		if (usedEdges.find(candidate) == usedEdges.end()) {
			return candidate;
		}
	}

	return -1;
}

std::vector<int> CuttingRailsImpl::SearchLoop(const int startIndex,
	const VertexToContinueIndexMap& vertexMap,
	const std::vector<ContinueEdgePtr>& toBeSearched,
	std::set<int>& usedEdges,
	const TopTools_IndexedMapOfShape& borderVertices,
	double squareTol) {
	std::vector<int> loop;
	loop.push_back(startIndex);
	usedEdges.insert(startIndex);
	// 获取起始连续边的端点
	TopoDS_Vertex startVertex = toBeSearched[startIndex]->GetStartVertex();
	gp_Pnt currentPoint = BRep_Tool::Pnt(toBeSearched[startIndex]->GetEndVertex());
	gp_Pnt startPoint = BRep_Tool::Pnt(startVertex);
	TopoDS_Shape startFace = toBeSearched[startIndex]->face;
	while (true) {
		// 查找下一条边
		int nextEdgeIndex = FindNextEdge(currentPoint, vertexMap, usedEdges, squareTol);

		if (nextEdgeIndex == -1) {
			// 没有找到下一条边，环不封闭
			loop.clear();
			break;
		}

		// 获取下一条边的端点
		TopoDS_Vertex nextVertex = toBeSearched[nextEdgeIndex]->GetEndVertex();
		TopoDS_Shape nextFace = toBeSearched[nextEdgeIndex]->face;
		gp_Pnt nextP = BRep_Tool::Pnt(nextVertex);

		// 确定下一个点
		currentPoint = nextP;

		// 添加到环中
		loop.push_back(nextEdgeIndex);
		usedEdges.insert(nextEdgeIndex);

		// 检查是否回到起点

		bool isSameFace = nextFace.IsSame(startFace);
		if (isSameFace && startVertex.IsSame(nextVertex)) {
			break;
		}

		if (!isSameFace) {
			double sqDist = currentPoint.SquareDistance(startPoint);
			if (sqDist <= squareTol && borderVertices.Contains(nextVertex)) {
				// 异面容差下共点且另一个为边界顶点
				break;
			}
		}
	}

	return loop;
}

std::vector<std::vector<CuttingRailsImpl::ContinueEdgePtr>> CuttingRailsImpl::PopLoop(const std::vector<int>& loop, const std::vector<ContinueEdgePtr>& totalContinues, double squareTol) {
	std::vector<std::vector<ContinueEdgePtr>> subLoops;
	std::vector<ContinueEdgePtr> subLoop;
	for (int i = 0, iLen = loop.size(); i < iLen; i++) {
		auto& curContinue = totalContinues[loop[i]];
		subLoop.push_back(curContinue);
		gp_Pnt end = curContinue->end;
		TopoDS_Shape nextFace = curContinue->face;
		TopoDS_Shape nextVertex = curContinue->GetEndVertex();
		for (int j = 0, jLen = subLoop.size() - 1; j < jLen; j++) {
			TopoDS_Shape startFace = subLoop[j]->face;
			bool isSameFace = nextFace.IsSame(startFace);
			bool canPop = false;
			if (isSameFace && subLoop[j]->GetStartVertex().IsSame(nextVertex)) {
				canPop = true;
			}
			if (!canPop && !isSameFace) {
				gp_Pnt start = subLoop[j]->start;
				double sqDist = start.SquareDistance(end);
				if (sqDist <= squareTol) {
					if (i != iLen - 1) {
						// 如果下一条更远，说明可以pop
						if (totalContinues[loop[i + 1]]->end.SquareDistance(start) > sqDist) {
							canPop = true;
						}
					}
					else {
						canPop = true;
					}

				}
				continue;
			}
			if (canPop) {
				std::vector<ContinueEdgePtr> loop(subLoop.begin() + j, subLoop.end());
				subLoops.push_back(loop);
				subLoop.erase(subLoop.begin() + j, subLoop.end());
				break;
			}
		}
	}
	if (!subLoop.empty()) {
		subLoops.push_back(subLoop);
	}
	return subLoops;
}

void CuttingRailsImpl::GetBorderVertex(const std::vector<ContinueEdgePtr>& continues, TopTools_IndexedMapOfShape& borderVertices) {
	borderVertices.Clear();
	for (const auto& continueEdge : continues) {
		if (continueEdge->degenerated) {
			continue;
		}
		// 获取起始和结束顶点
		TopoDS_Vertex startVertex = continueEdge->GetStartVertex();
		TopoDS_Vertex endVertex = continueEdge->GetEndVertex();
		// 添加到边界顶点映射中
		if (!borderVertices.Contains(startVertex)) {
			borderVertices.Add(startVertex);
		}
		else {
			borderVertices.RemoveKey(startVertex);
		}
		if (!borderVertices.Contains(endVertex)) {
			borderVertices.Add(endVertex);
		}
		else {
			borderVertices.RemoveKey(endVertex);
		}
	}
}

std::vector<CuttingRailsImpl::ContinueEdgePtr> CuttingRailsImpl::FindLoops(const std::vector<ContinueEdgePtr>& continues, std::vector<std::vector<ContinueEdgePtr>>& loops, const double tolerance) {
	// 将已成环的连续边pop成loop
	loops.clear();
	const double squareTol = tolerance * tolerance;
	std::vector<ContinueEdgePtr> toBeSearched;
	for (auto& curContinue : continues) {
		if (!curContinue->degenerated && (curContinue->GetStartVertex().IsSame(curContinue->GetEndVertex()))) {
			if (curContinue->len > 1.0e-6) {
				loops.push_back({ curContinue });
			}
		}
		else {
			if (!curContinue->degenerated) {
				toBeSearched.push_back(curContinue);
			}
		}
	}

	// 正常情况下不需要去重, 进入到搜环的连续边是唯一的, 端点处重合应该只会存在两条, 不需要转向选择
	TopTools_IndexedMapOfShape borderVertices;
	GetBorderVertex(toBeSearched, borderVertices);
	//auto vmap = showTopos(borderVertices);
	std::set<int> usedContinues;
	std::set<int> remainings;
	for (int i = 0, iLen = toBeSearched.size(); i < iLen; i++) {
		remainings.insert(i);
	}
	// 创建顶点到边的映射
	VertexToContinueIndexMap vertexMap = BuildVertexToContinueIndexMap(toBeSearched);

	// 遍历所有边，寻找可能的起始边
	for (int i = 0, iLen = toBeSearched.size(); i < iLen; i++) {
		const ContinueEdgePtr startContinue = toBeSearched[i];

		// 如果边已经被使用，跳过
		if (usedContinues.find(i) != usedContinues.end()) {
			continue;
		}

		// 从这条边开始搜索环
		std::vector<int> loop = SearchLoop(i, vertexMap, toBeSearched, usedContinues, borderVertices, squareTol);
		if (!loop.empty()) {
			std::vector<ContinueEdgePtr> edgeLoop;
			for (auto& id : loop) {
				edgeLoop.push_back(toBeSearched[id]);
				remainings.erase(id);
			}
			loops.push_back(edgeLoop);
		}
	}
	std::vector<ContinueEdgePtr> remainContinues;
	for (auto& id : remainings) {
		remainContinues.push_back(toBeSearched[id]);
	}
	return remainContinues;
}

TopoDS_Face CuttingRailsImpl::FixFaceWiresOrder(const TopoDS_Face& face) {
	// 创建新面
	BRep_Builder builder;
	TopoDS_Face newFace = TopoDS::Face(face.EmptyCopied());

	TopExp_Explorer wireExp(face, TopAbs_WIRE);

	while (wireExp.More()) {
		TopoDS_Wire innerWire = TopoDS::Wire(wireExp.Current());

		Handle(ShapeFix_Wire) fixInnerWire = new ShapeFix_Wire(innerWire, face, BRep_Tool::Tolerance(face));
		fixInnerWire->ModifyTopologyMode() = Standard_True;

		fixInnerWire->FixReorder();
		builder.Add(newFace, fixInnerWire->Wire());

		wireExp.Next();
	}

	return newFace;
}

TopoDS_Shape CuttingRailsImpl::FixFaceWiresOrder(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeShape& faceMapping) {
	TopoDS_Shell shell;
	BRep_Builder BB;
	BB.MakeShell(shell);

	// Clear the mapping first
	faceMapping.Clear();

	TopExp_Explorer faceExp(shape, TopAbs_FACE);

	while (faceExp.More()) {
		const TopoDS_Face& currentFace = TopoDS::Face(faceExp.Current());
		TopoDS_Face fixedFace = FixFaceWiresOrder(currentFace);

		// Store the mapping between original and fixed face
		faceMapping.Add(currentFace, fixedFace);

		BB.Add(shell, fixedFace);
		faceExp.Next();
	}

	return shell;
}

TopoDS_Shape CuttingRailsImpl::RemoveShortEdges(const TopoDS_Shape& originalShape, Standard_Real minLength) {
	// 用于重定向 (或移除) 拓扑实体的工具
	ShapeBuild_ReShape reshaper;

	// 遍历 originalShape 中所有 Edge
	for (TopExp_Explorer exp(originalShape, TopAbs_EDGE); exp.More(); exp.Next()) {
		const TopoDS_Edge& edge = TopoDS::Edge(exp.Current());
		if (BRep_Tool::Degenerated(edge)) {
			continue;
		}

		// 拿到几何曲线并计算长度
		BRepAdaptor_Curve adaptorCurve(edge);

		GProp_GProps props;
		BRepGProp::LinearProperties(edge, props);
		double length = props.Mass();

		// 如果长度小于指定阈值，就移除
		if (length < minLength) {
			reshaper.Remove(edge);

		}
	}

	// 应用 ReShape，生成新的 Shape
	TopoDS_Shape newShape = reshaper.Apply(originalShape);

	return newShape;
}

// 平面返回normal, 拉伸面返回 dir
bool CuttingRailsImpl::IsExtrudeFace(const TopoDS_Face& face, gp_Vec& dir1, gp_Vec& dir2, double theAngularTolerance) {
	// 获取面上的几何曲线
	BRepAdaptor_Surface surface(face);
	if (surface.GetType() == GeomAbs_Plane) {
		dir1.SetXYZ(surface.Plane().XAxis().Direction().XYZ());
		dir2.SetXYZ(surface.Plane().YAxis().Direction().XYZ());
		return true;
	}
	else if (surface.GetType() == GeomAbs_Cylinder) {
		dir1.SetXYZ(surface.Cylinder().Axis().Direction().XYZ());
		dir2.SetXYZ(surface.Cylinder().Axis().Direction().XYZ());
		return true;
	}
	else if (surface.GetType() == GeomAbs_SurfaceOfExtrusion) {
		dir1.SetXYZ(surface.Direction().XYZ());
		dir2.SetXYZ(surface.Direction().XYZ());
		return true;
	}
	else if (surface.GetType() == GeomAbs_BSplineSurface) {
		Standard_Integer udegree = surface.UDegree();
		Standard_Integer vdegree = surface.VDegree();
		Standard_Real u0 = surface.FirstUParameter();
		Standard_Real u1 = surface.LastUParameter();
		Standard_Real v0 = surface.FirstVParameter();
		Standard_Real v1 = surface.LastVParameter();
		if (udegree == 1 && vdegree == 1) {
			// 通过 D1() 来获取 U 方向的一阶导数向量
			Standard_Real um = 0.5 * (u0 + u1);
			Standard_Real vm = 0.5 * (v0 + v1);

			gp_Pnt pt;
			gp_Vec dU, dV;
			surface.D1(um, vm, pt, dU, dV);

			dU.Normalize();
			dV.Normalize();
			dir1.SetXYZ(dV.XYZ());  // 把 dU 视为拉伸方向
			dir2.SetXYZ(dU.XYZ());  // 把 dU 视为拉伸方向

			//if (!dV.IsEqual(gp_Vec(0, 0, 0), 1.0e-7, 1.0e-11)) {
			//	dV.Normalize();
			//	dir1.SetXYZ(dV.XYZ());  // 把 dU 视为拉伸方向
			//}
			return true;
		}
		else if (udegree == 1) {
			// 通过 D1() 来获取 U 方向的一阶导数向量
			Standard_Real um = 0.5 * (u0 + u1);
			Standard_Real vm = 0.5 * (v0 + v1);
			Handle(Geom_BSplineSurface) bsplineSurface = surface.BSpline();
			gp_Pnt pt;
			gp_Vec dU, dV;
			surface.D1(um, vm, pt, dU, dV);

			if (!dU.IsEqual(gp_Vec(0, 0, 0), 1.0e-7, 1.0e-11)) {
				dU.Normalize();
				std::vector<gp_Pnt> vControls;
				for (int iLen = surface.NbVPoles(), i = 1; i <= iLen; i++) {
					vControls.push_back(bsplineSurface->Pole(1, i));
				}
				gp_Vec normal;
				if (IsPlanerCurve(vControls, normal) > 0 && dU.IsParallel(normal, theAngularTolerance)) {
					dir1.SetXYZ(dU.XYZ());  // 把 dU 视为拉伸方向
					dir2.SetXYZ(dU.XYZ());  // 把 dU 视为拉伸方向
				}
				else {
					return false;
				}

			}
			return true;
		}
		else if (vdegree == 1) {
			Standard_Real um = 0.5 * (u0 + u1);
			Standard_Real vm = 0.5 * (v0 + v1);
			Handle(Geom_BSplineSurface) bsplineSurface = surface.BSpline();
			gp_Pnt pt;
			gp_Vec dU, dV;
			surface.D1(um, vm, pt, dU, dV);

			if (!dV.IsEqual(gp_Vec(0, 0, 0), 1.0e-7, 1.0e-11)) {
				dV.Normalize();

				std::vector<gp_Pnt> uControls;
				for (int j = 1, jLen = surface.NbUPoles(); j <= jLen; j++) {
					uControls.push_back(bsplineSurface->Pole(j, 1));
				}
				gp_Vec normal;
				if (IsPlanerCurve(uControls, normal) > 0 && dV.IsParallel(normal, theAngularTolerance)) {
					dir1.SetXYZ(dV.XYZ());  // 把 dV 视为拉伸方向
					dir2.SetXYZ(dV.XYZ());  // 把 dV 视为拉伸方向
				}
				else {
					return false;
				}
			}
			return true;
		}
	}
	else if (surface.GetType() == GeomAbs_SurfaceOfRevolution) {
		auto curve = surface.BasisCurve();
		gp_Pnt p;
		gp_Vec dirA;
		gp_Vec dirB;
		gp_Vec dirAB;
		curve->D1(curve->FirstParameter(), p, dirA);
		curve->D1(curve->LastParameter(), p, dirB);
		curve->D1((curve->FirstParameter() + curve->LastParameter()) / 2.0, p, dirAB);
		if (dirA.IsParallel(surface.AxeOfRevolution().Direction(), 0.001) &&
			dirB.IsParallel(surface.AxeOfRevolution().Direction(), 0.001) &&
			dirAB.IsParallel(surface.AxeOfRevolution().Direction(), 0.001)) {
			dir1.SetXYZ(surface.AxeOfRevolution().Direction().XYZ());  // 把 dV 视为拉伸方向
			dir2.SetXYZ(surface.AxeOfRevolution().Direction().XYZ());  // 把 dV 视为拉伸方向
			return true;
		}
	}
	return false;
}

double CuttingRailsImpl::CalculateAreaExteriorOnly(const TopoDS_Face& face) {
	// 获取外环
	TopoDS_Wire outerWire = BRepTools::OuterWire(face);

	// 创建只有外环的新面
	BRepBuilderAPI_MakeFace makeFace(BRep_Tool::Surface(face), outerWire);
	TopoDS_Face newFace = makeFace.Face();

	// 计算面积
	GProp_GProps gprops;
	BRepGProp::SurfaceProperties(newFace, gprops);
	return fabs(gprops.Mass());
}

TopoDS_Shape CuttingRailsImpl::SortFaceByArea(const TopoDS_Shape& shape, std::vector<double>& areas) {
	areas.clear();
	TopExp_Explorer explorer;
	std::vector<std::pair<TopoDS_Face, double>> sortedFaces;

	// 遍历所有face并存储面积
	for (explorer.Init(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
		TopoDS_Face currentFace = TopoDS::Face(explorer.Current());
		GProp_GProps gprops;
		BRepGProp::SurfaceProperties(currentFace, gprops);
		double area = fabs(gprops.Mass());
		sortedFaces.push_back(std::make_pair(currentFace, area));
	}

	// 按面积降序排序
	sort(sortedFaces.begin(), sortedFaces.end(),
		[](const std::pair<TopoDS_Face, double>& a, const std::pair<TopoDS_Face, double>& b) {
			return a.second > b.second;
		});

	BRep_Builder BB;
	TopoDS_Compound compound;
	BB.MakeCompound(compound);

	for (const auto& face : sortedFaces)
	{
		BB.Add(compound, face.first);
		areas.push_back(face.second);
	}
	return compound;
}

TopoDS_Shape CuttingRailsImpl::SortShellByArea(const TopoDS_Shape& shape) {
	TopExp_Explorer explorer;
	std::vector<std::pair<TopoDS_Shell, double>> sortedShells;

	// 遍历所有face并存储面积
	for (explorer.Init(shape, TopAbs_SHELL); explorer.More(); explorer.Next()) {
		TopoDS_Shell currentShell = TopoDS::Shell(explorer.Current());
		GProp_GProps gprops;
		BRepGProp::SurfaceProperties(currentShell, gprops);
		double area = fabs(gprops.Mass());
		sortedShells.push_back(std::make_pair(currentShell, area));
	}

	// 按面积降序排序
	sort(sortedShells.begin(), sortedShells.end(),
		[](const std::pair<TopoDS_Shell, double>& a, const std::pair<TopoDS_Shell, double>& b) {
			return a.second > b.second;
		});

	BRep_Builder BB;
	TopoDS_Compound compound;
	BB.MakeCompound(compound);

	for (const auto& shell : sortedShells)
	{
		BB.Add(compound, shell.first);
	}
	return compound;
}

void CuttingRailsImpl::ComputeFaceNormals(
	const TopoDS_Shape& shape,
	gp_Dir& normalStart,    // at (umin,vmin)
	gp_Dir& normalMiddle,   // at (u50,v50)
	gp_Dir& normalEnd)      // at (u100,v100)
{
	// 1. 将 shape 转成 Face
	TopoDS_Face face = TopoDS::Face(shape);

	// 2. 获取 Face 的底层几何曲面
	Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
	if (surface.IsNull())
	{
		// 若无法拿到 Surface，则返回默认值
		normalStart = gp::DZ();
		normalMiddle = gp::DZ();
		normalEnd = gp::DZ();
		return;
	}

	// 3. 获取 face 在 UV 范围的最小值和最大值
	Standard_Real uMin, uMax, vMin, vMax;
	BRepTools::UVBounds(face, uMin, uMax, vMin, vMax);

	// 计算 u 和 v 方向的中点参数值
	Standard_Real uLen = uMax - uMin;
	Standard_Real vLen = vMax - vMin;
	Standard_Real u50 = uMin + 0.50 * uLen;
	Standard_Real v50 = vMin + 0.50 * vLen;

	// 声明一个小函数计算 (u, v) 处的法向向量 (未归一化)
	auto computeOneNormal = [&](Standard_Real u, Standard_Real v) -> gp_Vec
		{
			gp_Pnt point;
			gp_Vec du, dv;
			surface->D1(u, v, point, du, dv);  // 计算一阶偏导
			return du.Crossed(dv);             // 叉乘得到法向
		};

	// 4. 分别计算三个点处的法向向量
	gp_Vec vecStart = computeOneNormal(uMin, vMin);    // 起点
	gp_Vec vecMiddle = computeOneNormal(u50, v50);     // 中点
	gp_Vec vecEnd = computeOneNormal(uMax, vMax);      // 终点

	// 5. 如果 face 的 Orientation 为 REVERSED，则需要反转法向
	if (face.Orientation() == TopAbs_REVERSED)
	{
		vecStart.Reverse();
		vecMiddle.Reverse();
		vecEnd.Reverse();
	}

	// 6. 归一化并写入引用参数
	vecStart.Normalize();
	vecMiddle.Normalize();
	vecEnd.Normalize();

	normalStart = gp_Dir(vecStart);
	normalMiddle = gp_Dir(vecMiddle);
	normalEnd = gp_Dir(vecEnd);
}

TopoDS_Shape CuttingRailsImpl::FixFaceOritation(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfos, TopTools_IndexedMapOfShape& reversedFaces)
{

	TopTools_IndexedMapOfShape visited;

	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {

		const TopoDS_Face benchmarkFace = TopoDS::Face(exp.Current());
		// 这里使用 queue 实现“真正的”BFS
		std::queue<TopoDS_Face> queue;
		queue.push(benchmarkFace);

		if (visited.Contains(benchmarkFace)) {
			continue;
		}

		// BFS 主循环
		while (!queue.empty())
		{
			// 取队列头
			TopoDS_Face curFace = TopoDS::Face(queue.front());
			queue.pop();

			// 若已经访问过，则跳过
			if (visited.Contains(curFace))
			{
				continue;
			}
			// 标记为已访问
			visited.Add(curFace);

			// 遍历 coedgeInfos 来找与 curFace 相邻的面
			for (auto& edge : coedgeInfos)
			{
				// 每项 edge.first 至少包含两条 coedge(面A, 面B)，edge.second 表示它们是否相反
				if (edge.first.size() != 2) {
					continue;
				}
				auto faceA = TopoDS::Face(edge.first[0]->face);
				auto faceB = TopoDS::Face(edge.first[1]->face);

				// 如果 faceA 就是当前面 curFace
				if (faceA.IsEqual(curFace))
				{
					bool faceAisReversed = reversedFaces.Contains(faceA);

					// edge.second == true => faceA,faceB 方向相反 => 需要翻转其中一个
					// 若 faceA 没翻转且方向相反 => 那 faceB 就得翻转
					// 若 faceA 已翻转且方向相反 => faceB 就无需翻转(二次翻转抵消)
					// => 这里可用异或: needFlipB = (faceAisReversed ^ edge.second)
					bool needFlipB = (faceAisReversed ^ edge.second);

					// 如果需要翻转 faceB，就加到 reversedFaces
					if (needFlipB)
					{
						reversedFaces.Add(faceB);
					}
					// 无论翻不翻，faceB 都要放进队列，以便后续 BFS
					queue.push(faceB);
				}
				// 如果 faceB 就是当前面
				else if (faceB.IsEqual(curFace))
				{
					bool faceBisReversed = reversedFaces.Contains(faceB);
					bool needFlipA = (faceBisReversed ^ edge.second);

					if (needFlipA)
					{
						reversedFaces.Add(faceA);
					}
					queue.push(faceA);
				}
			}
		}
	}

	TopoDS_Shell compound;
	BRep_Builder BB;
	BB.MakeShell(compound);

	TopExp_Explorer expFace(shape, TopAbs_FACE);
	for (; expFace.More(); expFace.Next()) {
		TopoDS_Face face = TopoDS::Face(expFace.Current());
		if (reversedFaces.Contains(face)) {
			face.Reverse();
		}
		BB.Add(compound, face);
	}
	return compound;
}

bool CuttingRailsImpl::ClassifyInOutGroups(const TopoDS_Shape& groups, TopoDS_Shape& outerBoundFaces, TopoDS_Shape& innerBoundFaces)
{
	BRep_Builder builder;
	TopoDS_Compound outers;
	TopoDS_Compound inners;
	builder.MakeCompound(outers);
	builder.MakeCompound(inners);

	int outerCount = 0;
	int innerCount = 0;

	for (TopExp_Explorer exp(groups, TopAbs_SHELL); exp.More(); exp.Next()) {
		TopoDS_Solid solid;
		builder.MakeSolid(solid);
		builder.Add(solid, exp.Current());
		BRepClass3d_SolidClassifier SC(solid);
		SC.PerformInfinitePoint(Precision::Confusion());
		if (SC.State() == TopAbs_OUT) {
			builder.Add(outers, exp.Current());
			outerCount++;
		}
		else {
			builder.Add(inners, exp.Current());
			innerCount++;
		}
	}

	TopoDS_Compound result;
	builder.MakeCompound(result);

	if (outerCount > 0 && innerCount > 0) {
		outerBoundFaces = outers;
		innerBoundFaces = inners;
	}
	else if (outerCount > 0) {
		outerBoundFaces = outers;
	}
	else if (innerCount > 0) {
		outerBoundFaces = inners;
	}
	else {
		return false;
	}
	
	return true;
}

bool HasCommonInMap(const TopTools_IndexedMapOfShape& map1, const TopTools_IndexedMapOfShape& map2) {
	for (int i = 1; i <= map1.Extent(); i++) {
		if (map2.Contains(map1.FindKey(i))) {
			return true;
		}
	}
	return false;
}

void CuttingRailsImpl::ClassifyInOutGroups(
	const TopoDS_Compound& groups,
	TopoDS_Shape& outerBoundFaces, TopoDS_Shape& innerBoundFaces, const bool hasExtrudeDir, const gp_Vec& extrudeDir, const SingleBorderType& singleBorderType)
{
	BRep_Builder builder;
	int borderCompCount = NbChildren(groups);
	TopoDS_Compound compInner, compOuter;
	builder.MakeCompound(compOuter);
	builder.MakeCompound(compInner);
	if (borderCompCount == 0) {
		outerBoundFaces = compOuter;
		innerBoundFaces = compInner;
	}
	else if (borderCompCount == 1) {
		if (singleBorderType == SingleBorderType::Inner) {
			outerBoundFaces = compOuter;
			innerBoundFaces = groups;
		}
		else {
			outerBoundFaces = groups;
			innerBoundFaces = compInner;
		}
	}
	else if (hasExtrudeDir) {
		// 采用投票的方式进行计算
		for (TopoDS_Iterator expShell(groups); expShell.More(); expShell.Next()) {
			auto& shell = expShell.Value();
			int status = ShapeInOutStatusByMesh(shell, extrudeDir);
			if (status == TopAbs_OUT) {
				builder.Add(compOuter, shell);
			}
			else {
				builder.Add(compInner, shell);
			}
		}
		if (NbChildren(compOuter) != 0 && NbChildren(compInner) != 0) {
			outerBoundFaces = compOuter;
			innerBoundFaces = compInner;
		}
	}
	// 面积比较兜底
	if (outerBoundFaces.IsNull() && innerBoundFaces.IsNull()) {
		TopoDS_Compound compOuter2, compInner2;
		builder.MakeCompound(compOuter2);
		builder.MakeCompound(compInner2);
		auto& biggestShell = GetBiggestGroup(groups);
		builder.Add(compOuter2, biggestShell);
		outerBoundFaces = compOuter2;
		TopExp_Explorer expInner(groups, TopAbs_SHELL);
		while (expInner.More()) {
			if (expInner.Current().IsNotEqual(biggestShell)) {
				builder.Add(compInner2, expInner.Current());
			}
			expInner.Next();
		}
		innerBoundFaces = compInner2;
	}

	return;
}

bool CuttingRailsImpl::IsShellConnect(const TopoDS_Shell& shell1, const TopoDS_Shell& shell2, const Face2CoedgesMap& face2Coedges) {
	for (TopoDS_Iterator exp(shell1); exp.More(); exp.Next()) {
		auto& face = exp.Value();

		auto& edges = face2Coedges.FindFromKey(face);
		for (const auto& edge : *edges)
		{
			if (edge.second && HasFaceInShell(TopoDS::Face(edge.second->face), shell2)) {
				return true;
			}
		}
	}
	return false;
}

void CuttingRailsImpl::ClassifyFuzzyFacesByInOutGroups(TopoDS_Shape& outerBoundFaces, TopoDS_Shape& innerBoundFaces,
	const TopoDS_Compound& fuzzyUnRecognized, const Face2CoedgesMap& face2Coedges)
{
	FaceIndexMap borderFacesMap;
	int index = 1;
	for (TopExp_Explorer exp(outerBoundFaces, TopAbs_SHELL); exp.More(); exp.Next()) {
		for (TopExp_Explorer expFace(exp.Current(), TopAbs_FACE); expFace.More(); expFace.Next()) {
			borderFacesMap.Add(expFace.Current(), index);
		}
		index++;
	}

	int outerBorderCount = index - 1;
	if (!innerBoundFaces.IsNull()) {
		for (TopExp_Explorer exp(innerBoundFaces, TopAbs_SHELL); exp.More(); exp.Next()) {
			for (TopExp_Explorer expFace(exp.Current(), TopAbs_FACE); expFace.More(); expFace.Next()) {
				borderFacesMap.Add(expFace.Current(), index);
			}
			index++;
		}
	}

	BRep_Builder builder;
	TopoDS_Compound outerFuzzy;
	TopoDS_Compound innerFuzzy;
	builder.MakeCompound(outerFuzzy);
	builder.MakeCompound(innerFuzzy);
	std::vector<TopoDS_Shell> outerAdded;
	std::vector<TopoDS_Shell> innerAdded;
	std::vector<bool> canFunzzyOuterAdd;
	std::vector<bool> canFunzzyInnerAdd;
	for (TopExp_Explorer expFuzzy(fuzzyUnRecognized, TopAbs_SHELL); expFuzzy.More(); expFuzzy.Next()) {
		auto& curFuzzy = expFuzzy.Current();
		std::set<int> connected;
		for (TopExp_Explorer expFuzzyFace(curFuzzy, TopAbs_FACE); expFuzzyFace.More(); expFuzzyFace.Next()) {
			const auto& edges = face2Coedges.FindFromKey(expFuzzyFace.Current());
			for (const auto& edge : *edges)
			{
				if (edge.second && borderFacesMap.Contains(edge.second->face)) {
					connected.insert(borderFacesMap.FindFromKey(edge.second->face));
				}
			}
		}

		if (connected.size() < 1) {
			continue;
		}

		int outerCount = 0;
		int innerCount = 0;
		for (auto connectIndex : connected) {
			connectIndex > outerBorderCount ? innerCount++ : outerCount++;
		}

		if (innerCount == 0) {
			outerAdded.push_back(TopoDS::Shell(curFuzzy));
			canFunzzyOuterAdd.push_back(true);
		}
		else if (outerCount == 0) {
			innerAdded.push_back(TopoDS::Shell(curFuzzy));
			canFunzzyInnerAdd.push_back(true);
		}
	}

	// 遍历outerAdded, 判断inner是否存在相邻, 如果存在则认为属于过渡面

	for (int i = 0, iLen = outerAdded.size(); i < iLen; i++) {
		auto& fuzzyOuter = outerAdded[i];
		for (int j = 0, jLen = innerAdded.size(); j < jLen; j++) {
			if (canFunzzyInnerAdd[j] && IsShellConnect(fuzzyOuter, innerAdded[j], face2Coedges)) {
				canFunzzyOuterAdd[i] = false;
				canFunzzyInnerAdd[j] = false;
				break;
			}
		}
	}

	for (int i = 0, iLen = outerAdded.size(); i < iLen; i++) {
		if (canFunzzyOuterAdd[i]) {
			builder.Add(outerBoundFaces, outerAdded[i]);
		}
	}
	
	for (int i = 0, iLen = innerAdded.size(); i < iLen; i++) {
		if (canFunzzyInnerAdd[i]) {
			builder.Add(innerBoundFaces, innerAdded[i]);
		}
	}
}

double CuttingRailsImpl::GetMinEdgeLength(const CoedgeInfo& coedgeInfos) {
	double minLength = 1000000.0;
	for (size_t i = 0; i < coedgeInfos.size(); i++) {
		auto& ce = coedgeInfos[i].first;
		double totalLength = 0.0;
		auto& edges = ce[0]->edges;
		for (auto& edge : edges)
		{
			{
				GProp_GProps gprops;
				BRepGProp::LinearProperties(edge, gprops);
				Standard_Real length = gprops.Mass();
				totalLength += length;
			}

		}

		if (minLength > totalLength) {
			minLength = totalLength;
		}
	}
	return minLength;

}

TopoDS_Compound CuttingRailsImpl::CreateCompoundFromFacePairs(
	const std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>>& facePairs)
{
	//====================================================
	// 第一步：收集所有出现过的 Face，并放入 shapeMap
	//====================================================
	TopTools_IndexedMapOfShape shapeMap;
	// shapeMap 是 1-based 索引容器

	for (const auto& pairFaces : facePairs)
	{
		shapeMap.Add(pairFaces.first);
		shapeMap.Add(pairFaces.second);
	}
	// shapeMap 现在包含了所有出现过的 Face（自动去重）

	//====================================================
	// 第二步：构造邻接表 adjacency
	//====================================================
	const Standard_Integer nbFaces = shapeMap.Extent(); // Face 的总数

	// 如果没有任何 face，则返回一个空的 compound
	if (nbFaces == 0)
	{
		TopoDS_Compound empty;
		BRep_Builder builder;
		builder.MakeCompound(empty);
		return empty;
	}

	// TColStd_Array1OfListOfInteger 是一个 1-based 的 array，
	// 下标从 1 到 nbFaces
	TColStd_Array1OfListOfInteger adjacency(1, nbFaces);

	// 根据输入的 pair，建立双向邻接
	// 注意：i1, i2 都是 1-based 的下标
	for (const auto& pairFaces : facePairs)
	{
		Standard_Integer i1 = shapeMap.FindIndex(pairFaces.first);
		Standard_Integer i2 = shapeMap.FindIndex(pairFaces.second);

		// 防御式检查：查找不到就会返回0，表示没在shapeMap里
		if (i1 == 0 || i2 == 0 || i1 == i2)
			continue;

		adjacency(i1).Append(i2);
		adjacency(i2).Append(i1);
	}

	//====================================================
	// 第三步：使用 BFS 找出所有连通分量
	//====================================================
	TColStd_PackedMapOfInteger visited; // 标记访问过的 Face 索引(1-based)

	// 每个连通分量存储它所包含的 face 索引
	std::vector<std::vector<Standard_Integer>> components;

	// 遍历所有 1-based face 索引
	for (Standard_Integer idx = 1; idx <= nbFaces; ++idx)
	{
		if (!visited.Contains(idx))
		{
			// 发现一个尚未访问的 Face 索引，将其所在的连通分量收集起来
			std::vector<Standard_Integer> currentComponent;

			std::queue<Standard_Integer> bfsQueue;
			bfsQueue.push(idx);
			visited.Add(idx);

			// 标准 BFS
			while (!bfsQueue.empty())
			{
				Standard_Integer curr = bfsQueue.front();
				bfsQueue.pop();
				currentComponent.push_back(curr);

				// 遍历与 curr 相邻的所有 face 索引
				// 注意要使用 TColStd_ListIteratorOfListOfInteger
				TColStd_ListIteratorOfListOfInteger it(adjacency(curr));
				for (; it.More(); it.Next())
				{
					Standard_Integer neigh = it.Value();
					if (!visited.Contains(neigh))
					{
						visited.Add(neigh);
						bfsQueue.push(neigh);
					}
				}
			}

			// 当前连通分量 BFS 完成
			components.push_back(currentComponent);
		}
	}

	//====================================================
	// 第四步：对每个连通分量创建一个 Shell，并加入顶层 Compound
	//====================================================
	BRep_Builder builder;
	TopoDS_Compound topCompound;
	builder.MakeCompound(topCompound);

	for (auto& comp : components)
	{
		TopoDS_Shell shell;
		builder.MakeShell(shell);

		// 将该连通分量内所有 Face 加入 Shell
		for (Standard_Integer faceIndex : comp)
		{
			// 通过索引获取原始 Face
			const TopoDS_Face& face = TopoDS::Face(shapeMap.FindKey(faceIndex));
			builder.Add(shell, face);
		}

		// 将 Shell 加入顶层 Compound
		builder.Add(topCompound, shell);
	}

	return topCompound;
}

TopoDS_Shape CuttingRailsImpl::MakeGroups(const CoedgeInfo& coedges, const TopTools_IndexedMapOfShape& inneroutFaces) {
	std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>> facePairs;
	for (size_t i = 0; i < coedges.size(); i++)
	{
		if (coedges[i].first.size() == 2) {
			auto face0 = coedges[i].first[0]->face;
			auto face1 = coedges[i].first[1]->face;

			if (inneroutFaces.Contains(face0) && inneroutFaces.Contains(face1)) {
				facePairs.push_back({ face0, face1 });
			}
		}
	}
	TopoDS_Shape compoundShells = CreateCompoundFromFacePairs(facePairs);
	return compoundShells;
}

TopoDS_Shape CuttingRailsImpl::RemoveInnerLoopFaces(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfos, double thin)
{
	if (thin > Precision::Confusion()) {
		return shape;
	}

	BRep_Builder builder;
	TopoDS_Shell newShell = TopoDS::Shell(shape.EmptyCopied());
	for (TopoDS_Iterator exp(shape); exp.More(); exp.Next()) {
		auto& face = TopoDS::Face(exp.Value());
		if (NbChildren(face) == 1) {
			builder.Add(newShell, face);
			continue;
		}
		TopoDS_Face newFace = TopoDS::Face(face.EmptyCopied());
		TopoDS_Wire outerWire = BRepTools::OuterWire(face);

		for (TopoDS_Iterator expWire(face); expWire.More(); expWire.Next()) {
			auto& wire = TopoDS::Wire(expWire.Value());
			if (wire.IsSame(outerWire)) {
				builder.Add(newFace, wire);
				continue;
			}
			bool isMatched = false;
			for (TopoDS_Iterator expEdge(wire); expEdge.More(); expEdge.Next()) {
				auto& edge = TopoDS::Edge(expEdge.Value());
				for (auto& coedges : coedgeInfos) {
					auto& matched = coedges.first;
					if (matched.size() == 1) {
						continue;
					}

					for (auto& continues : matched) {
						if (continues->face.IsSame(face) && !continues->bOuterLoop) {
							isMatched = std::any_of(continues->edges.begin(), continues->edges.end(), [&](const TopoDS_Edge& e) {
								return e.IsSame(edge);
								});
						}
						if (isMatched) {
							break;
						}
					}

					if (isMatched) {
						break;
					}
				}

				if (isMatched) {
					break;
				}
			}

			if (isMatched) {
				builder.Add(newFace, wire);
			}
		}
		
		builder.Add(newShell, newFace);
	}
	return newShell;
}

bool CuttingRailsImpl::GetFaceInnerPt(const TopoDS_Face& face, double& u, double& v)
{
	BRepClass_FaceExplorer explore(face);
	BRepClass_FClassifier classifier;
	double tol = BRep_Tool::Tolerance(face);
	double umin, umax, vmin, vmax;

	TopoDS_Wire outer;
	TopExp_Explorer expWire(face, TopAbs_WIRE);
	if (NbChildren(face) > 1) {
		outer = BRepTools::OuterWire(face);
	}
	else {
		outer = TopoDS::Wire(expWire.Current());
	}
	int gpMaxSteps = 32;
	
	BRepTools::UVBounds(face, umin, umax, vmin, vmax);
	double uLen = umax - umin;
	double vLen = vmax - vmin;

	double innerX = 0;
	double innerY = 0;

	bool found = false;
	
	for (int steps = 1; steps <= gpMaxSteps && !found; steps *= 2) {
		double uHalfstep = uLen / (steps * 2.0);
		double vHalfstep = vLen / (steps * 2.0);

		for (int i = 0; i < steps && !found; ++i) {
			innerX = umin + uHalfstep * (1 + 2 * i);

			for (int j = 0; j < steps && !found; ++j) {
				innerY = vmin + vHalfstep * (1 + 2 * j);
				classifier.Perform(explore, gp_Pnt2d(innerX, innerY), tol);
				found = classifier.State() == TopAbs_IN;
			}
		}
	}

	if (found) {
		u = innerX;
		v = innerY;
		return true;
	}

	// Find a point near the contour geomLoop
	double rLen = std::min(uLen / (gpMaxSteps * 4.0), vLen / (gpMaxSteps * 4.0));
	
	for (TopExp_Explorer exp(outer, TopAbs_EDGE); exp.More(); exp.Next()) {
		TopoDS_Edge edge = TopoDS::Edge(exp.Current());
		if (BRep_Tool::Degenerated(edge)) {
			continue;
		}
		Standard_Real firstParam, lastParam;
		Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(edge, face, firstParam, lastParam);
		if (pcurve.IsNull()) {
			continue;
		}
		Standard_Real midParam = 0.5 * (firstParam + lastParam);

		std::vector<double> refParams = { midParam, edge.Orientation() == TopAbs_FORWARD ? firstParam : lastParam };
		for (auto& param : refParams) {
			// Gradually approaching to the reference point to find the inner point
			gp_Vec2d tangent;
			gp_Pnt2d uvPt;
			pcurve->D1(param, uvPt, tangent);
			if (tangent.SquareMagnitude() <= gp::Resolution()) {
				continue;
			}
			tangent.Normalize();
			double tangentX = tangent.X();
			double tangentY = tangent.Y();
			for (int k = 0; k < gpMaxSteps; ++k) {
				// Looking in multiple directions
				double curRlen = rLen / std::pow(2, k);
				gp_Vec2d normal(-tangentY * curRlen, tangentX * curRlen);
				gp_Vec2d quaderDir = normal.Added(tangent).Normalized().Multiplied(curRlen);
				std::vector<gp_Pnt2d> judgePoints = {
				uvPt.Translated(normal),
				uvPt.Translated(normal.Reversed()),
				uvPt.Translated(quaderDir),
				uvPt.Translated(quaderDir.Reversed())};

				for (auto& judgePoint : judgePoints) {
					classifier.Perform(explore, judgePoint, tol);
					if (classifier.State() == TopAbs_IN) {
						u = judgePoint.X();
						v = judgePoint.Y();
						return true;
					}
				}
			}
		}

	}
	return false;
}

void CuttingRailsImpl::GetConnectFaceSeamEdges(const std::vector<std::pair<ContinueEdgePtr, bool>>& coedges, const CoedgeInfo& coedgeInfos, std::vector<ContinueEdgePtr>& seamEdges)
{
	auto& face = TopoDS::Face(coedges[0].first->face);
	if (NbChildren(face) != 1) {
		return;
	}
	std::set<ContinueEdgePtr> otherContinues;
	for (auto& infos : coedgeInfos) {
		auto& coedges = infos.first;
		if (coedges.size() == 2) {
			if (coedges[0]->face.IsSame(face) && !coedges[0]->degenerated) {
				otherContinues.insert(coedges[0]);
			}
			else if (coedges[1]->face.IsSame(face) && !coedges[1]->degenerated) {
				otherContinues.insert(coedges[1]);
			}
		}
	}
	
	for (auto& cur : coedges) {
		otherContinues.erase(cur.first);
	}

	if (otherContinues.size() == 0) {
		return;
	}

	int startIndex = seamEdges.size();
	std::vector<std::pair<gp_Lin, bool>> directions;
	for (auto& cur : coedges) {
		bool sameDir = cur.second;
		auto& curContinue = cur.first;

		auto& toBeFindVertex = sameDir ? TopExp::FirstVertex(curContinue->edges.front(), true) : TopExp::LastVertex(curContinue->edges.back(), true);
		for (auto& single : otherContinues) {
			TopoDS_Vertex V1 = sameDir ? TopExp::LastVertex(single->edges.back(), true) : TopExp::FirstVertex(single->edges.front(), true);

			// 检查是否为起点
			if (V1.IsSame(toBeFindVertex)) {
				// 当single为线性时, 认为是seam edge
				gp_Lin line(single->start, gp_Vec(single->start, single->end));
				double t = ElCLib::Parameter(line, single->mid);
				double squareDist = ElCLib::Value(t, line).SquareDistance(single->mid);
				if (squareDist < 1.0e-12) {
					seamEdges.push_back(single);
					directions.push_back({ line, sameDir });
				}
				otherContinues.erase(single);
				break;
			}
		}
	}
	
	int endIndex = seamEdges.size();
	// 缝合边可能不完整, 需要找回
	for (int i = startIndex, j = 0; i < endIndex; i++, j++) {

		std::vector<ContinueEdgePtr> toBeMerge;
		toBeMerge.push_back(seamEdges[i]);
		int count = 0;
		int maxCount = otherContinues.size();
		auto& line = directions[j].first;
		bool sameDir = directions[j].second;
		ContinueEdgePtr cur = seamEdges[i];
		while (cur && count++ < maxCount) {
			auto& toBeFindVertex = sameDir ? TopExp::FirstVertex(cur->edges.front(), true) : TopExp::LastVertex(cur->edges.back(), true);
			cur = nullptr;
			for (auto& single : otherContinues) {
				TopoDS_Vertex V1 = sameDir ? TopExp::LastVertex(single->edges.back(), true) : TopExp::FirstVertex(single->edges.front(), true);

				// 检查是否为起点
				if (V1.IsSame(toBeFindVertex)) {
					// 当single为线性时, 认为是seam edge
					double t1 = ElCLib::Parameter(line, single->mid);
					double squareDist1 = ElCLib::Value(t1, line).SquareDistance(single->mid);
					gp_Pnt refPt = sameDir ? single->start : single->end;
					double t2 = ElCLib::Parameter(line, refPt);
					double squareDist2 = ElCLib::Value(t2, line).SquareDistance(refPt);
					if (squareDist1 < 1.0e-12 && squareDist2 < 1.0e-12) {
						cur = single;
						toBeMerge.push_back(single);
					}
					otherContinues.erase(single);
					break;
				}
			}
		}

		if (toBeMerge.size() > 1) {
			double totalLen = 0;
			std::vector<TopoDS_Edge> es;
			for (auto& single : toBeMerge) {
				totalLen += single->len;
				for (auto& e : single->edges) {
					es.push_back(e);
				}
			}
			seamEdges[i] = std::make_shared<ContinueEdge>(es, totalLen);
		}
	}
	return;
}


void CuttingRailsImpl::GetConnectFaceSeamEdges(const TopoDS_Face& face, const CoedgeInfo& coedgeInfos, std::vector<ContinueEdgePtr>& seamEdges)
{
	if (NbChildren(face) == 1) {
		return;
	}
	std::vector<std::pair<ContinueEdgePtr, bool>> faceInnerConnectContinues;
	NCollection_IndexedDataMap<TopoDS_Shape, std::vector<std::pair<ContinueEdgePtr, bool>>> FaceCountMap;

	for (auto& infos : coedgeInfos) {
		auto& coedges = infos.first;
		if (coedges.size() == 2) {
			auto& coedge1 = coedges[0];
			auto& coedge2 = coedges[1];
			if (!coedge1->bOuterLoop && coedge2->bOuterLoop && coedge1->face.IsSame(face)) {
				faceInnerConnectContinues.push_back({ coedge2, infos.second });
				if (FaceCountMap.Contains(coedge2->face)) {
					FaceCountMap.ChangeFromKey(coedge2->face).push_back({ coedge2, infos.second });
				}
				else {
					FaceCountMap.Add(coedge2->face, { { coedge2, infos.second } });
				}
			}
			else if (!coedge2->bOuterLoop && coedge1->bOuterLoop && coedge2->face.IsSame(face)) {
				faceInnerConnectContinues.push_back({ coedge1, infos.second });
				if (FaceCountMap.Contains(coedge1->face)) {
					FaceCountMap.ChangeFromKey(coedge1->face).push_back({ coedge1, infos.second });
				}
				else {
					FaceCountMap.Add(coedge1->face, { { coedge1, infos.second } });
				}
			}
		}
	}

	for (int i = 1; i <= FaceCountMap.Extent(); i++) {
		GetConnectFaceSeamEdges(FaceCountMap.FindFromIndex(i), coedgeInfos, seamEdges);
	}
}

int CuttingRailsImpl::SolidInOutStatus(const TopoDS_Solid& solid, bool isSingle, const gp_Vec& dir)
{
	if (isSingle) {
		// 单壁处理
		return ShapeInOutStatusByMesh(solid, dir);
	}

	BRepClass3d_SolidClassifier SC(solid);
	SC.PerformInfinitePoint(Precision::Confusion());

	return SC.State();
}

bool CuttingRailsImpl::GetSortedShapeExtrudeDir(const TopoDS_Shape& shape, const std::vector<double>& sortedFaceAreas, gp_Vec& dir)
{
	bool isExtrude = GetExtrudeDirection(shape, dir);
	if (isExtrude) {
		return true;
	}
	// std::cout << "没有找到拉伸方向, 尝试使用近似拉伸计算" << std::endl;
	bool find = GetNearExtrudeDirection(shape, sortedFaceAreas, dir);
	if (!find) {
		find = FindDirectionBySurfaceNormalCross(shape, dir);
	}

	return find;
}

int CuttingRailsImpl::ShapeInOutStatusByMesh(const TopoDS_Shape& shape, const gp_Vec& dir)
{
	double xmin, xmax, ymin, ymax, zmin, zmax;
	gp_Cylinder cylinder = GetAproxCylinder(shape, dir, xmin, ymin, zmin, xmax, ymax, zmax);
	return ShapeInOutStatusByMesh(shape, dir, cylinder, xmin, ymin, zmin, xmax, ymax, zmax);
}

gp_Cylinder CuttingRailsImpl::GetAproxCylinder(const TopoDS_Shape& shape, const gp_Vec& dir, double& xmin, double& ymin, double& zmin, double& xmax, double& ymax, double& zmax) {
	// 生成一个圆柱面
	TopTools_IndexedMapOfShape vertices;
	TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
	Bnd_Box box;
	for (int i = 1; i <= vertices.Extent(); i++) {
		box.Add(BRep_Tool::Pnt(TopoDS::Vertex(vertices(i))));
	}

	box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

	gp_Pnt center((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2);
	gp_Ax2 cylAx2(center, dir);
	double radius = std::max(xmax - xmin, std::max(ymax - ymin, zmax - zmin));
	return gp_Cylinder(cylAx2, radius);
}

int CuttingRailsImpl::ShapeInOutStatusByMesh(const TopoDS_Shape& shape, const gp_Vec& dir, const gp_Cylinder& cylinder, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
{
	// 对shape做离散, 如果shape本就存在离散, 注意orientation
	{
		const double THE_LINEAR_DEFLECTION = 5.e-3;
		const double THE_ANGULAR_DEFLECTION = 0.5;
		const bool   THE_IS_PARALLEL = false;

		BRepMesh_FastDiscret::Parameters aParameters;
		const gp_XYZ          aMin(xmin, ymin, zmin);
		const gp_XYZ          aMax(xmax, ymax, zmax);
		const double          aDeflection = std::max((aMax - aMin).Modulus() * THE_LINEAR_DEFLECTION, THE_LINEAR_DEFLECTION);
		aParameters.Deflection = aDeflection;
		aParameters.MinSize = aDeflection * 2;
		aParameters.Angle = THE_ANGULAR_DEFLECTION;
		aParameters.InParallel = THE_IS_PARALLEL;
		aParameters.ControlSurfaceDeflection = true;
		// aParameters.EnableControlSurfaceDeflectionAllSurfaces = true;
		BRepMesh_IncrementalMesh mesher(shape, aParameters);
		mesher.Perform();
	}

	// 筛选出靠近圆柱侧面的面片, 将面片中点投影到圆柱面上, 判断面片中点法向和投影点法向是否一致，如果一致，则是朝外, 如果不一致, 则是朝内
	// 统计：多少三角片朝内，多少三角片朝外
	int triCountOuter = 0;
	int triCountInner = 0;

	// 阈值：15°
	double angleThreshold = cos(M_PI / 12);
	for (TopExp_Explorer fexp(shape, TopAbs_FACE); fexp.More(); fexp.Next())
	{
		TopoDS_Face face = TopoDS::Face(fexp.Current());
		
		// Face 的离散网格
		TopLoc_Location loc;
		Handle(Poly_Triangulation) polyTri = BRep_Tool::Triangulation(face, loc);
		if (polyTri.IsNull()) {
			continue; // 没有网格
		}
		// face.Orientation() 会影响网格法向的正反
		bool isReversed = (face.Orientation() == TopAbs_REVERSED);

		int nbTriangle = polyTri->NbTriangles();
		// 遍历所有三角形
		for (Standard_Integer iTri = 1; iTri <= nbTriangle; ++iTri)
		{
			// tri 是一个三角形，有3个顶点索引
			Poly_Triangle tri = polyTri->Triangle(iTri);
			Standard_Integer n1, n2, n3;
			tri.Get(n1, n2, n3);
			gp_Pnt p1 = polyTri->Node(n1);
			gp_Pnt p2 = polyTri->Node(n2);
			gp_Pnt p3 = polyTri->Node(n3);

			// 计算三角形法向
			//      normal = (p2 - p1) x (p3 - p1)
			gp_Vec v12(p1, p2);
			gp_Vec v13(p1, p3);
			gp_Vec triNormal = v12.Crossed(v13);
			// 如果 face 是 REVERSED，需要反转
			if (isReversed) {
				triNormal.Reverse();
			}

			// Check 长度是否足够大，以免零长度出错
			Standard_Real normMag = triNormal.Magnitude();
			if (normMag < 1e-10) {
				// degenerate triangle
				continue;
			}
			triNormal.Divide(normMag); // 归一化

			// 4.2) 判断 triNormal 和 shapeDir 的夹角
			double dotValue = triNormal.Dot(dir);
			if (std::abs(dotValue) > angleThreshold) {
				// 小于 60° => 贴近“轴向”，不在圆柱侧面，不处理
				continue;
			}
			// >= 60° => 说明这个三角面片“偏离轴向明显”，接近“圆柱侧面” 

			// 4.3) 取三角形中点
			gp_Pnt ptMid(
				(p1.X() + p2.X() + p3.X()) / 3.0,
				(p1.Y() + p2.Y() + p3.Y()) / 3.0,
				(p1.Z() + p2.Z() + p3.Z()) / 3.0
			);

			// 计算 Cylinder 的参数
			gp_Pnt projOnCyl;
			gp_Vec dU, dV;
			double u, v;
			ElSLib::Parameters(cylinder, ptMid, u, v);
			ElSLib::D1(u, v, cylinder, projOnCyl, dU, dV);

			// 投影点上计算圆柱面的法向
			gp_Vec cylNormal = dU.Crossed(dV);
			cylNormal.Normalize();

			// 4.6) 比较 triNormal 与 cylNormal 的点积
			double dotVal = triNormal.Dot(cylNormal);
			// dotVal > 0 => 同向 => 可视作“外”
			// dotVal < 0 => 反向 => 可视作“内”
			if (dotVal >= 0.0) {
				++triCountOuter;
			}
			else {
				++triCountInner;
			}
		} // end for iTri
	} // end for fexp
	if (triCountOuter > triCountInner)
		return TopAbs_OUT;  // maybe "大部分朝外"
	else
		return TopAbs_IN; // maybe "大部分朝内"
}

TopoDS_Shape CuttingRailsImpl::GetBiggestGroup(const TopoDS_Shape& groups) {
	TopExp_Explorer exp(groups, TopAbs_SHELL);
	double maxArea = -1.0;
	TopoDS_Shape ret;
	for (; exp.More(); exp.Next())
	{
		GProp_GProps props;
		BRepGProp::SurfaceProperties(exp.Current(), props);
		Standard_Real area = props.Mass();

		if (area > maxArea) {
			maxArea = area;
			ret = exp.Current();
		}
	}

	return ret;
}

TopoDS_Shape CuttingRailsImpl::ShowEdges(const TopTools_MapOfShape& edges) {
	BRep_Builder builder;
	TopoDS_Compound compound;
	builder.MakeCompound(compound);

	// 遍历 TopTools_MapOfShape
	for (TopTools_MapIteratorOfMapOfShape it(edges); it.More(); it.Next())
	{
		const TopoDS_Shape& aShape = it.Key();
		// 将这个 shape 加入到 compound 中
		builder.Add(compound, aShape);
	}

	return compound;
}

TopoDS_Shape CuttingRailsImpl::ShowEdges(const TopTools_IndexedMapOfShape& edges) {
	TopoDS_Compound compound;
	BRep_Builder builder;
	builder.MakeCompound(compound);

	// Assuming you have a TopTools_IndexedMapOfShape called 'map'
	for (Standard_Integer i = 1; i <= edges.Extent(); i++) {
		const TopoDS_Shape& shape = edges.FindKey(i);
		builder.Add(compound, shape);
	}

	return compound;
}

TopoDS_Shape CuttingRailsImpl::RemoveInnerLoopFaces(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfos, TopTools_IndexedMapOfShape& innerLoopFaces, double thin) {
	;
	for (size_t i = 0; i < coedgeInfos.size(); i++)
	{
		if (coedgeInfos[i].first.size() == 2) {
			auto& continuesEdgesA = coedgeInfos[i].first[0];
			auto& continuesEdgesB = coedgeInfos[i].first[1];

			if (!continuesEdgesA->bOuterLoop) {
				if (NbChildren(continuesEdgesA->face) == 2) {
					if (CDistance(TopoDS::Face(continuesEdgesA->face)) > 1.1 * thin) {
						innerLoopFaces.Add(continuesEdgesB->face);
					}
				}
				else {
					innerLoopFaces.Add(continuesEdgesB->face);
				}
			}

			if (!continuesEdgesB->bOuterLoop) {
				if (NbChildren(continuesEdgesB->face) == 2) {
					if (CDistance(TopoDS::Face(continuesEdgesB->face)) > 1.1 * thin) {
						innerLoopFaces.Add(continuesEdgesA->face);
					}
				}
				else {
					innerLoopFaces.Add(continuesEdgesA->face);
				}
			}
		}
	}

	// Create a new compound to store the result
	BRep_Builder builder;
	TopoDS_Compound result;
	builder.MakeCompound(result);

	// Traverse all faces in the shape
	TopExp_Explorer explorer(shape, TopAbs_FACE);
	for (; explorer.More(); explorer.Next()) {
		const TopoDS_Face& face = TopoDS::Face(explorer.Current());

		// If the face is not in innerLoopFaces, add it to the result
		if (!innerLoopFaces.Contains(face)) {
			builder.Add(result, face);
		}
		else {
			continue;
		}
	}

	return result;
}

Standard_Real CuttingRailsImpl::ComputeThicknessByRay(
	const TopoDS_Shape& shape, const TopoDS_Face& skipFace, const CoedgeInfo& coedgeInfoVector)
{
	std::vector<ContinueEdgePtr> seamEdges;
	GetConnectFaceSeamEdges(skipFace, coedgeInfoVector, seamEdges);

	double seamEdgeAvgLength = 0;

	int seamEdgesCount = seamEdges.size();
	if (seamEdgesCount > 0) {
		for (auto& seam : seamEdges) {
			seamEdgeAvgLength += seam->len;
		}
		seamEdgeAvgLength /= seamEdgesCount;
	}

	if (seamEdgeAvgLength > 0.5 && seamEdgeAvgLength < 2) {
		return seamEdgeAvgLength;
	}

	gp_Pnt pt;
	gp_Vec du;
	gp_Vec dv;
	double aU, aV;
	bool find = false;
	gp_Vec normal;
	if (GetFaceInnerPt(skipFace, aU, aV)) {
		BRep_Tool::Surface(skipFace)->D1(aU, aV, pt, du, dv);
		normal = du.Crossed(dv);
		if (normal.SquareMagnitude() > gp::Resolution()) {
			find = true;
		}
	}

	if(!find) {
		BRepClass3d_SolidExplorer aSE(skipFace);
		double aParam = 0.5;
		aSE.FindAPointInTheFace(skipFace, pt, aU, aV, aParam, du, dv);
		normal = du.Crossed(dv);
	}
	
	if (skipFace.Orientation() == TopAbs_REVERSED) {
		normal.Reverse();
	}

	gp_Lin rayLine(pt, normal);
	//auto rayEdge = BRepBuilderAPI_MakeEdge(rayLine, -20, 20).Edge();

	Standard_Real minDist = std::numeric_limits<Standard_Real>::max();
	
	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next())
	{
		TopoDS_Face face = TopoDS::Face(exp.Current());

		// 可根据需求跳过与同一面的相交检测
		if (face.IsSame(skipFace))
			continue;

		// 构造 IntCurvesFace_ShapeIntersector
		IntCurvesFace_ShapeIntersector intersector;
		// face, line, tolerance
		intersector.Load(face, Precision::Confusion());
		intersector.Perform(rayLine, -std::numeric_limits<Standard_Real>::max(), std::numeric_limits<Standard_Real>::max());

		if (intersector.IsDone() && intersector.NbPnt() > 0)
		{
			// 获取所有交点，找距离最小者
			for (Standard_Integer i = 1; i <= intersector.NbPnt(); ++i)
			{
				Standard_Real dist = std::abs(intersector.WParameter(i));
				if (dist < minDist)
				{
					minDist = dist;
				}
			}
		}
	}

	if (minDist == std::numeric_limits<Standard_Real>::max())
	{
		// 表示未找到有效交点
		return -1.0;
	}
	else
	{
		return minDist;
	}
}

bool CuttingRailsImpl::IsLeftRightOnSameSide(const EdgeVectorPtr& coedges, const ContinueEdgePtr& coedge, double tol) {

	bool bStartSame = true;
	bool bEndSame = true;
	for (const auto& current : *coedges)
	{
		if (current.first == coedge) {
			continue;
		}

		if (current.first->degenerated) {
			continue;
		}

		bool bStartMath = current.first->end.IsEqual(coedge->start, tol);
		bool bEndMath = current.first->start.IsEqual(coedge->end, tol);
		if (bStartMath) {
			gp_Vec V1 = current.first->dirAtEnd();
			gp_Vec V2 = coedge->dirAtStart();

			gp_Vec normal;
			if (current.first->NormalAt(coedge->uvAtStart(), normal)) {
				double val = normal.DotCross(V1, V2);

				if (val < -1.0e-6) {
					bStartSame = false;
				}
			}
			else {
				return true;
			}


		}

		if (bEndMath) {
			gp_Vec V1 = coedge->dirAtEnd();
			gp_Vec V2 = current.first->dirAtStart();

			gp_Vec normal;
			if (current.first->NormalAt(coedge->uvAtEnd(), normal)) {
				double val = normal.DotCross(V1, V2);

				if (val < -1.0e-6) {
					bEndSame = false;
				}
			}
			else {
				return true;
			}
		}
	}

	return bStartSame == bEndSame;
}

bool CuttingRailsImpl::IsSameClassification(const std::pair<ContinueEdgePtr, ContinueEdgePtr>& edge, const Face2CoedgesMap& faceCoedgesMap, double thinDist, double tol) {
	bool smooth = edge.first->edgetType == SMOOTH && edge.first->len > 2.0 * thinDist /*&&
		IsLeftRightOnSameSide(faceCoedgesMap.FindFromKey(edge.first->face), edge.first, tol) &&
		IsLeftRightOnSameSide(faceCoedgesMap.FindFromKey(edge.second->face), edge.second, tol)*/;
	//IsFacesOnDifferentSide(edge, faceCoedgesMap, tol) true是同一类，flase其实是凹， 所以不需要这个条件了
	return smooth || (edge.first->edgetType == AO && edge.first->len > 2.0 * thinDist);
}

bool CuttingRailsImpl::IsSameClassificationWithoutLength(const std::pair<ContinueEdgePtr, ContinueEdgePtr>& edge, const Face2CoedgesMap& faceCoedgesMap, double tol) {
	bool smooth = edge.first->edgetType == SMOOTH /*&&
		IsLeftRightOnSameSide(faceCoedgesMap.FindFromKey(edge.first->face), edge.first, tol) &&
		IsLeftRightOnSameSide(faceCoedgesMap.FindFromKey(edge.second->face), edge.second, tol)*/;
	//IsFacesOnDifferentSide(edge, faceCoedgesMap, tol) true是同一类，flase其实是凹， 所以不需要这个条件了
	return smooth || (edge.first->edgetType == AO);
}

bool CuttingRailsImpl::IsFacesOnDifferentSide(const std::pair<ContinueEdgePtr, ContinueEdgePtr>& edge, const Face2CoedgesMap& map, double tol) {

	const auto& pt1 = edge.first->start;

	bool bStartMath = false;
	bool bEndMath = false;
	gp_Vec V1, V2;

	const auto& edges = map.FindFromKey(edge.first->face);
	for (const auto& current : *edges)
	{
		bEndMath = current.first->end.IsEqual(pt1, tol);
		if (bEndMath) {
			V1 = current.first->dirAtEnd();
			break;
		}
	}

	const auto& edges2 = map.FindFromKey(edge.second->face);
	for (const auto& current : *edges2)
	{
		bStartMath = current.first->start.IsEqual(pt1, tol);
		if (bStartMath) {
			V2 = current.first->dirAtStart();
			break;
		}
	}
	if (bStartMath && bEndMath) {
		if (V1.Dot(V2) > 0.0) {
			return true;
		}
		else {
			return false;
		}
	}

	return false;
}

bool CuttingRailsImpl::JudgeTwoWileFaceBeCuttingFace(const TopoDS_Shape& face, const EdgeVectorPtr& faceCoedges, double thinDist) {

	double f, l;

	std::vector<std::pair<ContinueEdgePtr, ContinueEdgePtr>> outerEdges;
	std::vector<std::pair<ContinueEdgePtr, ContinueEdgePtr>> innerEdges;

	for (const auto& edge : *faceCoedges) {
		if (edge.first->bOuterLoop) {
			outerEdges.push_back(edge);
		}
		else {
			innerEdges.push_back(edge);
		}
	}

	std::sort(innerEdges.begin(), innerEdges.end(),
		[](const std::pair<ContinueEdgePtr, ContinueEdgePtr>& lhs,
			const std::pair<ContinueEdgePtr, ContinueEdgePtr>& rhs)
		{
			return lhs.first->len > rhs.first->len;
		});

	bool bMatch = false;

	for (const auto& innerEdge : innerEdges)
	{
		if (!innerEdge.second) {
			continue;
		}
		double uInnerAdj, vInnerAdj;
		const auto& innerAdjEdge = innerEdge.second->edges;
		const auto& midEdge = innerAdjEdge[innerAdjEdge.size() / 2];
		const auto& innerAdjFace = TopoDS::Face(innerEdge.second->face);
		Handle(Geom2d_Curve) innerAdjcurve2d = BRep_Tool::CurveOnSurface(midEdge, innerAdjFace, f, l);
		const Handle(Geom2d_TrimmedCurve)& trimedCurve = new Geom2d_TrimmedCurve(innerAdjcurve2d, f, l);
		const Handle(Geom_CurveOnSurface)& curveOnsurface = new Geom_CurveOnSurface(trimedCurve, BRep_Tool::Surface(innerAdjFace));
		//param = edge.Orientation() == face.Orientation() ? l : f;
		double param = (f + l) / 2;
		gp_Pnt2d pt;
		innerAdjcurve2d->D0(param, pt);
		uInnerAdj = pt.X();
		vInnerAdj = pt.Y();

		gp_Pnt midPt;
		gp_Vec midV;
		curveOnsurface->D1(param, midPt, midV);
		Handle(Geom_Plane) plane = new Geom_Plane(midPt, midV);

		// 使用 BRepBuilderAPI_MakeVertex 将点转换为拓扑顶点
		BRepBuilderAPI_MakeVertex mkVertex(midPt);
		TopoDS_Vertex midPtVertex = mkVertex.Vertex();

		double A = midV.X();
		double B = midV.Y();
		double C = midV.Z();
		double D = midV.X() * midPt.X() + midV.Y() * midPt.Y() + midV.Z() * midPt.Z();
		D = -D;

		double u = 0;
		double v = 0;
		double w = 0;
		double outAdjMinu = 0;
		double outAdjMinv = 0;
		double dist2 = 0;
		double minDist2 = 1000000.0;
		TopoDS_Face outerMinAdjFace;

		std::pair<ContinueEdgePtr, ContinueEdgePtr> targetEdge;
		for (const auto& outerEdge : outerEdges)
		{
			//真实求交
			if (!outerEdge.second) {
				continue;
			}
			bMatch = true;
			const auto& outerAdjFace = TopoDS::Face(outerEdge.second->face);

			auto edges = outerEdge.second->edges;
			for (const auto& edge : edges) {
				//bool bSameParam = BRep_Tool::SameParameter(TopoDS::Edge(exp.Current()));
				Handle(Geom2d_Curve) curve2d = BRep_Tool::CurveOnSurface(TopoDS::Edge(edge), outerAdjFace, f, l);
				const Handle(Geom2d_TrimmedCurve) trimedCurve = new Geom2d_TrimmedCurve(curve2d, f, l);
				const Handle(Geom_CurveOnSurface)& curveOnsurface = new Geom_CurveOnSurface(trimedCurve, BRep_Tool::Surface(outerAdjFace));

				GeomAPI_IntCS intCS(curveOnsurface, plane);
				if (intCS.IsDone()) {
					int nbPoints = intCS.NbPoints();
					for (int i = 1; i <= nbPoints; ++i) {
						gp_Pnt pt = intCS.Point(i);
						BRepBuilderAPI_MakeVertex mkVertex(pt);
						TopoDS_Vertex ptVertex = mkVertex.Vertex();

						intCS.Parameters(i, u, v, w);
						dist2 = u * u + v * v;
						if (minDist2 > dist2) {
							minDist2 = dist2;
							gp_Pnt2d outPt;
							curve2d->D0(w, outPt);
							outAdjMinu = outPt.X();
							outAdjMinv = outPt.Y();
							outerMinAdjFace = outerAdjFace;
							targetEdge = outerEdge;
						}
					}
				}
			}
		}

		if (std::sqrt(minDist2) < 1.5 * thinDist) {

			if (innerEdge.first->edgetType == AO && targetEdge.first->edgetType == AO) {
				return false;
			}

			gp_Dir aDNInnerAdj;
			gp_Dir aDNOuterAdj;
			FaceNormal(innerAdjFace, uInnerAdj, vInnerAdj, aDNInnerAdj);
			FaceNormal(outerMinAdjFace, outAdjMinu, outAdjMinv, aDNOuterAdj);
			double angle = aDNInnerAdj.Angle(aDNOuterAdj);
			if (angle > 0.8 * M_PI) {
				return true;
			}
		}
		else {
			continue;
		}

	}

	if (!bMatch) {
		const auto& midEdge = innerEdges[0].first->edges[innerEdges[0].first->edges.size() / 2];
		Handle(Geom_Curve) innerCurve = BRep_Tool::Curve(TopoDS::Edge(midEdge), f, l);
		//param = edge.Orientation() == face.Orientation() ? l : f;
		double param = (f + l) / 2;
		gp_Pnt midPt;
		gp_Vec midV;
		innerCurve->D1(param, midPt, midV);
		Handle(Geom_Plane) plane = new Geom_Plane(midPt, midV);

		double u, v, w, dist2;
		double minDist2 = 10000.0;
		for (const auto& outerEdge : outerEdges)
		{


			auto edges = outerEdge.first->edges;
			for (const auto& edge : edges) {
				//bool bSameParam = BRep_Tool::SameParameter(TopoDS::Edge(exp.Current()));
				Handle(Geom_Curve) curve = BRep_Tool::Curve(TopoDS::Edge(edge), f, l);

				GeomAPI_IntCS intCS(curve, plane);
				if (intCS.IsDone()) {
					int nbPoints = intCS.NbPoints();
					for (int i = 1; i <= nbPoints; ++i) {
						gp_Pnt pt = intCS.Point(i);
						BRepBuilderAPI_MakeVertex mkVertex(pt);
						TopoDS_Vertex ptVertex = mkVertex.Vertex();

						intCS.Parameters(i, u, v, w);
						dist2 = u * u + v * v;
						if (minDist2 > dist2) {
							minDist2 = dist2;
						}
					}
				}
			}
		}

		if (std::sqrt(minDist2) < 1.1 * thinDist) {
			return true;
		}

		else {
			return false;
		}

	}

	//判断不出来
	return false;
}

void CuttingRailsImpl::GetContinuesEdges(const TopoDS_Shape& face, const CoedgeInfo& coedges, std::vector<std::pair<ContinueEdgePtr, ContinueEdgePtr>>& faceCoedges) {
	faceCoedges.clear();
	for (auto& coedge : coedges) {
		if (coedge.first.size() == 2) {
			auto ceA = coedge.first[0];
			auto ceB = coedge.first[1];
			if (ceA->face.IsEqual(face)) {
				faceCoedges.push_back({ ceA, ceB });
			}
			else if (ceB->face.IsEqual(face)) {
				faceCoedges.push_back({ ceB, ceA });
			}
		}
		else if (coedge.first.size() == 1) {
			auto ceA = coedge.first[0];
			if (ceA->face.IsEqual(face)) {
				faceCoedges.push_back({ ceA, nullptr });
			}
		}
	}
}

void CuttingRailsImpl::AddFacesToGroup(const TopoDS_Shape& shape, TopTools_IndexedMapOfShape& faceGroup) {
	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
		faceGroup.Add(exp.Current());
	}
}

void CuttingRailsImpl::AddAdjecentSmoothFaces(const TopoDS_Shape& face, const Face2CoedgesMap& faceCoedgesMap,
	TopTools_IndexedMapOfShape& faceGroup, double tol) {


	// 这里使用 queue 实现“真正的”BFS
	std::queue<TopoDS_Shape> queue;
	queue.push(face);

	// BFS 主循环
	while (!queue.empty())
	{
		// 取队列头
		TopoDS_Face curFace = TopoDS::Face(queue.front());
		queue.pop();

		// 若已经访问过，则跳过
		if (faceGroup.Contains(curFace))
		{
			continue;
		}
		// 标记为已访问
		faceGroup.Add(curFace);

		// 遍历 coedgeInfos 来找与 curFace 相邻的面
		const auto& twinCoedges = faceCoedgesMap.FindFromKey(curFace);
		for (auto& edge : *twinCoedges)
		{
			if (!edge.second) {
				continue;
			}

			// 如果 faceA 就是当前面 curFace
			if (edge.first->edgetType == SMOOTH &&
				/*IsLeftRightOnSameSide(faceCoedgesMap.FindFromKey(edge.first->face), edge.first, tol) &&
				IsLeftRightOnSameSide(faceCoedgesMap.FindFromKey(edge.second->face), edge.second, tol) &&*/
				IsFacesOnDifferentSide(edge, faceCoedgesMap, tol))
			{
				queue.push(edge.second->face);
			}
		}
	}
}

TopoDS_Compound CuttingRailsImpl::GetSmoothFaces(const TopoDS_Face& face, const Face2CoedgesMap& faceCoedgesMap, double tol) {
	TopTools_IndexedMapOfShape faceGroup;

	// 这里使用 queue 实现“真正的”BFS
	std::queue<TopoDS_Shape> queue;
	queue.push(face);

	// BFS 主循环
	while (!queue.empty())
	{
		// 取队列头
		TopoDS_Face curFace = TopoDS::Face(queue.front());
		queue.pop();

		// 若已经访问过，则跳过
		if (faceGroup.Contains(curFace))
		{
			continue;
		}
		// 标记为已访问
		faceGroup.Add(curFace);

		// 遍历 coedgeInfos 来找与 curFace 相邻的面
		const auto& twinCoedges = faceCoedgesMap.FindFromKey(curFace);
		for (auto& edge : *twinCoedges)
		{
			if (!edge.second) {
				continue;
			}

			// 如果 faceA 就是当前面 curFace
			if (edge.first->edgetType == SMOOTH &&
				/*IsLeftRightOnSameSide(faceCoedgesMap.FindFromKey(edge.first->face), edge.first, tol) &&
				IsLeftRightOnSameSide(faceCoedgesMap.FindFromKey(edge.second->face), edge.second, tol) &&*/
				IsFacesOnDifferentSide(edge, faceCoedgesMap, tol))
			{
				queue.push(edge.second->face);
			}
		}
	}

	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);
	for (size_t i = 1; i <= faceGroup.Extent(); i++)
	{
		BB.Add(comp, faceGroup.FindKey(i));
	}

	return comp;
}

/**
 * @brief 查找仅出现一次的连续区间（若有多个，仅返回第一个满足条件的区间）。
 *
 * @param arr      输入数组
 * @param startIdx 通过引用返回的区间起始位置
 * @param endIdx   通过引用返回的区间结束位置
 * @return true    找到只出现一次的连续区间
 * @return false   未找到只出现一次的连续区间
 */
bool CuttingRailsImpl::findSingleRun(const std::vector<int>& arr, int& startIdx, int& endIdx) {
	// 如果数组为空，直接返回 false
	if (arr.empty()) {
		return false;
	}

	// 1. 收集所有连续区间
	std::vector<RunInfo> runs;
	int runStart = 0;  // 连续段起始
	for (int i = 1; i < static_cast<int>(arr.size()); ++i) {
		if (arr[i] != arr[i - 1]) {
			// 前后数字不相等，说明上一个连续段结束
			runs.push_back({ arr[i - 1], runStart, i - 1 });
			// 下一段起始位置
			runStart = i;
		}
	}
	// 处理最后一个连续段
	runs.push_back({ arr.back(), runStart, static_cast<int>(arr.size()) - 1 });

	// 2. 统计每个数字出现了多少个连续“段”
	//    （如数字 1 出现了哪几段，数字 2 又出现了哪几段等等）
	std::unordered_map<int, int> runCount;
	for (auto& run : runs) {
		runCount[run.value]++;
	}

	// 3. 在所有连续区间中，找出只出现过一次的那段（只出现过一段的数字）
	//    如果有多个，只返回最先找到的一个
	for (auto& run : runs) {
		if (runCount[run.value] == 1) {
			// 将目标区间的起止下标通过引用返回
			startIdx = run.start;
			endIdx = run.end;
			return true;
		}
	}

	// 没有找到只出现一次的连续区间
	return false;
}

TopoDS_Shape CuttingRailsImpl::GetHolesOtherSideFaces(const TopoDS_Shape& holes, const TopoDS_Shape& smoothFaces, const Face2CoedgesMap& faceCoedgesMap, double tol) {
	TopTools_IndexedMapOfShape faces;
	TopTools_IndexedMapOfShape otherSizeFaces;

	for (TopExp_Explorer exp(holes, TopAbs_FACE); exp.More(); exp.Next()) {
		const TopoDS_Face& face = TopoDS::Face(exp.Current());
		faces.Add(face);
	}

	for (TopExp_Explorer exp(smoothFaces, TopAbs_FACE); exp.More(); exp.Next()) {
		const TopoDS_Face& face = TopoDS::Face(exp.Current());
		faces.Add(face);
	}

	for (TopExp_Explorer exp(holes, TopAbs_FACE); exp.More(); exp.Next()) {
		const TopoDS_Face& face = TopoDS::Face(exp.Current());
		const auto edges = faceCoedgesMap.FindFromKey(face);
		for (const auto& edge : *edges)
		{
			if (edge.second && !faces.Contains(edge.second->face)) {
				otherSizeFaces.Add(edge.second->face);
			}
		}
	}

	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);

	for (size_t i = 1; i <= otherSizeFaces.Extent(); i++)
	{
		BB.Add(comp, otherSizeFaces.FindKey(i));
	}

	return comp;
}

TopoDS_Face CuttingRailsImpl::GetOtherSideFace(const TopoDS_Shape& sommothFace, BRepClass3d_SolidExplorer& aSE, double& minDist) {
	TopExp_Explorer exp(sommothFace, TopAbs_FACE);
	TopoDS_Face face = TopoDS::Face(exp.Current());
	Standard_Real  aU = 0., aV = 0.;
	minDist = 10000.0;
	gp_Pnt aPoint;
	gp_Dir aDN;

	double aParam = 0.5;
	bool bFound = aSE.FindAPointInTheFace(face, aPoint, aU, aV, aParam);
	FaceNormal(face, aU, aV, aDN);
	gp_Lin aLin(aPoint, aDN);

	Handle(Geom_Line) geomLine = new Geom_Line(aLin);
	TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(geomLine, -20.0, 20.0);

	TopoDS_Face theOtherFace;
	for (aSE.InitShell(); aSE.MoreShell(); aSE.NextShell()) {
		if (aSE.RejectShell(aLin) == Standard_False) {
			for (aSE.InitFace(); aSE.MoreFace(); aSE.NextFace()) {
				if (aSE.RejectFace(aLin) == Standard_False) {
					TopoDS_Shape aLocalShape = aSE.CurrentFace();
					TopoDS_Face CurFace = TopoDS::Face(aLocalShape);

					if (CurFace.IsSame(face)) {
						continue;
					}
					IntCurvesFace_Intersector& Intersector3d = aSE.Intersector(CurFace);
					Intersector3d.Perform(aLin, -200.0, 200.0);

					if (Intersector3d.IsDone()) {
						// 取该面与射线所有交点中“距离 p3d 最近的”点
						for (Standard_Integer i = 1; i <= Intersector3d.NbPnt(); ++i) {
							gp_Pnt ipnt = Intersector3d.Pnt(i);
							Standard_Real dist = ipnt.Distance(aPoint);
							if (dist < minDist && dist > Precision::Confusion()) {
								// 更新当前最近值
								minDist = dist;
								theOtherFace = CurFace;
							}
						}
					}
				}
			}
		}
	}

	return theOtherFace;
}

bool CuttingRailsImpl::IsInnerOuterShell(const TopoDS_Shell& shell, const Face2CoedgesMap& faceCoedgesMap, double thinDist, double tol) {

	//存在两个以上的内环
	int wileSize = 0;
	for (TopoDS_Iterator exp(shell); exp.More(); exp.Next()) {
		const TopoDS_Face& face = TopoDS::Face(exp.Value());
		wileSize += NbChildren(face);
	}
	if (wileSize >= NbChildren(shell) + 2) {
		return true; // 多环
	}

	// 搜出来的跨多个面的内环
	//TopoDS_Shape holes = GetSmoothFacesInnerLoopsFace(shell, faceCoedgesMap, tol);
	//if (holes.NbChildren() >= 2) {
	//	return true;
	//}

	// 面积/周长
	// 计算面积/周长 与 ratio的关系
	GProp_GProps gprops;
	BRepGProp::SurfaceProperties(shell, gprops);
	double area = fabs(gprops.Mass());

	GProp_GProps props;
	BRepGProp::LinearProperties(shell, props);
	double length = props.Mass();
	// 计算周长
	for (TopExp_Explorer exp(shell, TopAbs_FACE); exp.More(); exp.Next()) {
		const auto& faceCoedges = faceCoedgesMap.FindFromKey(exp.Current());
		for (const auto& edge : *faceCoedges) {
			if (IsSameClassification(edge, faceCoedgesMap, thinDist, tol)) {
				length -= edge.first->len;
			}
		}
	}

	double ratio = area / length * 2;

	if (ratio > 2.0 * thinDist) {
		return true;
	}

	return false;
}

bool CuttingRailsImpl::HasFaceInShell(const TopoDS_Face& targetFace, const TopoDS_Shell& shell)
{
	for (TopExp_Explorer exp(shell, TopAbs_FACE); exp.More(); exp.Next())
	{
		if (exp.Current().IsSame(targetFace))
		{
			return true;
		}
	}
	return false;
}

TopoDS_Shell CuttingRailsImpl::GetRelationPatchs(const TopoDS_Shape& otherSideFace, const TopoDS_Shape& fixedShape) {

	const TopoDS_Face& face = TopoDS::Face(otherSideFace);
	for (TopoDS_Iterator expShell(fixedShape); expShell.More(); expShell.Next())
	{
		if (HasFaceInShell(face, TopoDS::Shell(expShell.Value())))
		{
			return TopoDS::Shell(expShell.Value());
		}
	}

	BRep_Builder BB;
	TopoDS_Shell shell;
	BB.MakeShell(shell);

	return shell;

}

void CuttingRailsImpl::GetAdjecentFaces(const TopoDS_Shape& shell, const Face2CoedgesMap& face2Coedges,
	TopTools_IndexedMapOfShape& adjecentFaces) {

	TopTools_IndexedMapOfShape facesSet;
	for (TopExp_Explorer exp(shell, TopAbs_FACE); exp.More(); exp.Next())
	{
		facesSet.Add(exp.Current());
	}

	for (int i = 1; i <= facesSet.Extent(); i++) {
		const TopoDS_Face& face = TopoDS::Face(facesSet.FindKey(i));
		const auto edges = face2Coedges.FindFromKey(face);
		for (const auto& edge : *edges)
		{
			if (edge.second && !facesSet.Contains(edge.second->face)) {
				adjecentFaces.Add(edge.second->face);
			}
		}
	}
}

void CuttingRailsImpl::GetAdjecentSmoothFaces(const TopoDS_Shape& shell, const Face2CoedgesMap& face2Coedges,
	TopTools_IndexedMapOfShape& adjecentSmoothFaces, double tol) {

	TopTools_IndexedMapOfShape facesSet;
	for (TopExp_Explorer exp(shell, TopAbs_FACE); exp.More(); exp.Next())
	{
		facesSet.Add(exp.Current());
	}

	for (int i = 1; i <= facesSet.Extent(); i++) {
		const TopoDS_Face& face = TopoDS::Face(facesSet.FindKey(i));
		const auto edges = face2Coedges.FindFromKey(face);
		for (const auto& edge : *edges)
		{
			if (edge.first->degenerated) {
				continue;
			}
			if (IsSameClassificationWithoutLength(edge, face2Coedges, tol)) {
				adjecentSmoothFaces.Add(edge.second->face);
			}
		}
	}
}

CuttingRailsImpl::AdjecentType CuttingRailsImpl::CheckFacesDistribution(const TopTools_IndexedMapOfShape& currentPatchAdjecentFaces,
	const TopoDS_Shell& outerBorderFaces,
	const TopoDS_Shell& innerBorderFaces)
{
	bool hasOuterFace = false;
	bool hasInnerFace = false;

	// 遍历 currentPatchAdjecentFaces 中的每个面
	for (Standard_Integer i = 1; i <= currentPatchAdjecentFaces.Extent(); i++)
	{
		const TopoDS_Face& face = TopoDS::Face(currentPatchAdjecentFaces.FindKey(i));

		// 检查是否在外边界面中
		if (hasOuterFace || HasFaceInShell(face, outerBorderFaces))
		{
			hasOuterFace = true;
		}

		// 检查是否在内边界面中
		if (hasInnerFace || HasFaceInShell(face, innerBorderFaces))
		{
			hasInnerFace = true;
		}

		// 如果已经找到了同时存在于内外边界的情况，可以提前返回
		if (hasOuterFace && hasInnerFace)
		{
			return BOTH;
		}
	}

	if (hasOuterFace && !hasInnerFace) {
		return OUTER;
	}
	else if (!hasOuterFace && hasInnerFace) {
		return INNER;
	}
	else if (hasOuterFace && hasInnerFace) {
		return BOTH;
	}
	else {
		return UNKNOWNAdjecentType;
	}
}

bool CuttingRailsImpl::HasCommonFace(const TopoDS_Shape& patchShell, const TopTools_IndexedMapOfShape& adjecentFaces) {
	for (TopExp_Explorer it(patchShell, TopAbs_FACE); it.More(); it.Next()) {

		if (adjecentFaces.Contains(it.Current())) {
			return true;
		}
	}
	return false;
}

void CuttingRailsImpl::GetSameClassificationFromFuzzyFaces(const TopoDS_Face& face, const TopTools_IndexedMapOfShape& fuzzyFaces, const Face2CoedgesMap& faceCoedgesMap, TopTools_IndexedMapOfShape& faceGroup, double thinDist, double tol) {

	// 这里使用 queue 实现“真正的”BFS
	std::queue<TopoDS_Shape> queue;
	queue.push(face);

	// BFS 主循环
	while (!queue.empty())
	{
		// 取队列头
		TopoDS_Face curFace = TopoDS::Face(queue.front());
		queue.pop();

		// 若已经访问过，则跳过
		if (faceGroup.Contains(curFace))
		{
			continue;
		}
		// 标记为已访问
		faceGroup.Add(curFace);

		// 遍历 coedgeInfos 来找与 curFace 相邻的面
		const auto& twinCoedges = faceCoedgesMap.FindFromKey(curFace);
		for (auto& edge : *twinCoedges)
		{
			if (!edge.second || edge.first->degenerated) {
				continue;
			}

			// 如果 faceA 就是当前面 curFace
			if (IsSameClassificationWithoutLength(edge, faceCoedgesMap, tol) && fuzzyFaces.Contains(edge.second->face))
				queue.push(edge.second->face);
		}
	}
}

TopoDS_Compound CuttingRailsImpl::MakeSameClassificationForFuzzyFaces(const TopTools_IndexedMapOfShape& fuzzyFaces, const Face2CoedgesMap& faceCoedgesMap, double thinDist, double tol) {

	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);

	TopTools_IndexedMapOfShape faceGroup;
	int indexStart = 1;
	for (int i = 1; i <= fuzzyFaces.Extent(); i++) {
		TopoDS_Face face = TopoDS::Face(fuzzyFaces.FindKey(i));

		if (faceGroup.Contains(face)) {
			continue;
		}

		GetSameClassificationFromFuzzyFaces(face, fuzzyFaces, faceCoedgesMap, faceGroup, thinDist, tol);
		TopoDS_Shell shell;
		BB.MakeShell(shell);
		for (size_t i = indexStart; i <= faceGroup.Extent(); i++)
		{
			BB.Add(shell, faceGroup.FindKey(i));
		}
		BB.Add(comp, shell);

		indexStart = faceGroup.Extent() + 1;
	}
	return comp;
}

TopoDS_Shell CuttingRailsImpl::GetPatchFromFace(const TopoDS_Compound& comp, const TopoDS_Face& face) {
	for (TopoDS_Iterator it(comp); it.More(); it.Next())
	{
		TopoDS_Shell shell = TopoDS::Shell(it.Value());
		bool has = HasFaceInShell(face, shell);

		if (has) {
			return shell;
		}
	}
	return TopoDS_Shell();
}

bool CuttingRailsImpl::IsExtrude(const TopTools_IndexedMapOfShape& faces, gp_Vec& dir) {
	bool firstFace = true;
	double theAngularTolerance = 0.01;
	std::vector<gp_Vec> candidates;
	for (int i = 1; i <= faces.Extent(); i++) {

		auto face = TopoDS::Face(faces(i));
		gp_Vec dir1, dir2;
		bool gotDirs = IsExtrudeFace(face, dir1, dir2, theAngularTolerance);
		if (!gotDirs) {
			return false;
		}

		if (firstFace) {
			candidates.push_back(dir1);
			candidates.push_back(dir2);
			firstFace = false;
		}
		else {
			// 后续面的处理：
			// 只保留与当前面方向匹配的候选
			std::vector<gp_Vec> newCands;
			for (auto& cand : candidates) {
				// 如果与dir1或dir2相同，则保留
				if (cand.IsParallel(dir1, theAngularTolerance) || cand.IsParallel(dir2, theAngularTolerance)) {
					newCands.push_back(cand);
				}
			}
			candidates.swap(newCands);
			if (candidates.empty()) {
				// 如果候选方向已经空了，说明没有公共方向
				return false;
			}
		}
	}

	// 遍历完所有面后，若候选集不为空，取第一个作为最终结果
	if (!candidates.empty()) {
		dir.SetXYZ(candidates.front().XYZ());
		return true;
	}

	// 没有可用候选
	return false;
}

void CuttingRailsImpl::GetInnerOuterShellByRatio(const TopoDS_Shape& fixedShape, const Face2CoedgesMap& faceCoedgesMap,
	TopTools_IndexedMapOfShape& innerouterFaces, TopTools_IndexedMapOfShape& cuttingFaces,
	TopTools_IndexedMapOfShape& fuzzyFaces,
	std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair,
	bool bExtrudeShape, double thinDist, double tol) {

	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);
	BRepClass3d_SolidExplorer aSE(fixedShape);
	int index = 0;
	double minDist;
	for (TopoDS_Iterator exp(fixedShape); exp.More(); exp.Next()) {
		const TopoDS_Shell& shell = TopoDS::Shell(exp.Value());
		index++;
		const auto& firstFace = GetFirstSubShape(shell);
		if (innerouterFaces.Contains(firstFace) || cuttingFaces.Contains(firstFace)) {
			continue;
		}

		bool bBorderShell = IsInnerOuterShell(shell, faceCoedgesMap, thinDist, tol);
		if (bBorderShell) {

			AddFacesToGroup(shell, innerouterFaces);
			BB.Add(comp, shell);
			TopoDS_Shape otherSideFace = GetOtherSideFace(shell, aSE, minDist);
			TopoDS_Shell otherSidePatch = GetRelationPatchs(otherSideFace, fixedShape);

			if (shell.IsEqual(otherSidePatch)) {
				continue;
			}
			
			if (!innerouterFaces.Contains(otherSideFace)) {
				AddFacesToGroup(otherSidePatch, innerouterFaces);
				BB.Add(comp, otherSidePatch);
			}

			mutuallyExclusivePair.push_back({ shell , otherSidePatch });
		}
		else {
			for (TopoDS_Iterator exp(shell); exp.More(); exp.Next()) {
				fuzzyFaces.Add(exp.Value());
			}
		}
	}

	// auto outerShape = showTopos(innerouterFaces);
	// auto cuttingShape = showTopos(cuttingFaces);
	// auto fuzzyShape = showTopos(fuzzyFaces);

	TopoDS_Compound fuzzyComp = MakeSameClassificationForFuzzyFaces(fuzzyFaces, faceCoedgesMap, thinDist, tol);
	//要确保从大到小
	index = 0;

	for (TopoDS_Iterator it(fuzzyComp); it.More(); it.Next())
	{
		TopoDS_Shell shell = TopoDS::Shell(it.Value());
		index++;
		TopTools_IndexedMapOfShape adjecentFaces;
		GetAdjecentFaces(shell, faceCoedgesMap, adjecentFaces);
		//互斥
		bool IsAdjecentFacesMutuallyExclusive = false;
		for (const auto pair : mutuallyExclusivePair)
		{
			//otherSideFace在互斥对的一端， 邻面在互斥对的另一端
			if (HasCommonFace(pair.first, adjecentFaces) && HasCommonFace(pair.second, adjecentFaces)) {
				IsAdjecentFacesMutuallyExclusive = true;
				break;
			}
		}
		if (IsAdjecentFacesMutuallyExclusive) {
			continue;
		}

		TopoDS_Face otherSideFace = GetOtherSideFace(shell, aSE, minDist);
		// if (minDist > 2 * thinDist) {
		// 	continue;
		// }

		bool bValid = false;
		TopoDS_Shell otherSidePatch;
		for (const auto pair : mutuallyExclusivePair)
		{
			//     |---|              ---|
			// ----|   |------           |-------
			//  
			// ----------------      ------------
			//小凸起: otherSideFace在互斥对的一端， 邻面在互斥对的另一端
			bool match1 = HasFaceInShell(otherSideFace, pair.first) && HasCommonFace(pair.second, adjecentFaces);
			bool match = match1 || (HasFaceInShell(otherSideFace, pair.second) && HasCommonFace(pair.first, adjecentFaces));
			if (match) {
				bValid = true;
				if (match1) {
					otherSidePatch = pair.first;
				}
				else {
					otherSidePatch = pair.second;
				}
				break;
			}
		}
		if (bValid) {

			AddFacesToGroup(shell, innerouterFaces);
			BB.Add(comp, shell);

			mutuallyExclusivePair.push_back({ shell , otherSidePatch });
			continue;
		}

		TopTools_IndexedMapOfShape adjecentSmoothFaces;
		GetAdjecentSmoothFaces(shell, faceCoedgesMap, adjecentSmoothFaces, tol);

		//和borderFaces光滑过度，则为壁面
		for (size_t i = 1; i <= adjecentSmoothFaces.Extent(); i++)
		{
			if (innerouterFaces.Contains(adjecentSmoothFaces.FindKey(i))) {
				TopoDS_Shell smoothPatch = GetPatchFromFace(comp, TopoDS::Face(adjecentSmoothFaces.FindKey(i)));

				AddFacesToGroup(shell, innerouterFaces);
				BB.Add(comp, shell);

				for (const auto& pair : mutuallyExclusivePair)
				{
					if (pair.first.IsEqual(smoothPatch)) {
						mutuallyExclusivePair.push_back({ shell , pair.second });
					}

					if (pair.second.IsEqual(smoothPatch)) {
						mutuallyExclusivePair.push_back({ shell , pair.first });

					}
				}
				break;
			}
		}
	}

	return;
}

bool CuttingRailsImpl::IsOuterFace(const TopoDS_Face& face, BRepClass3d_SolidExplorer& aSE) {
	// iteratively try up to 10 probing points from each face
	Standard_Real  aU = 0., aV = 0.;
	gp_Pnt aPoint;
	gp_Dir aDN;
	const Standard_Integer NB_MAX_POINTS_PER_FACE = 10;
	int intersectNum = 0;
	for (Standard_Integer itry = 0; itry < NB_MAX_POINTS_PER_FACE; itry++)
	{
		double aParam = 0.5;
		bool bFound = aSE.FindAPointInTheFace(face, aPoint, aU, aV, aParam);
		if (!bFound || !FaceNormal(face, aU, aV, aDN))
			continue;
		gp_Lin aLin(aPoint, aDN);

		Handle(Geom_Line) geomLine = new Geom_Line(aLin);
		TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(geomLine, 0, 1000);

		intersectNum = 0;
		for (aSE.InitShell(); aSE.MoreShell(); aSE.NextShell()) {
			if (aSE.RejectShell(aLin) == Standard_False) {
				for (aSE.InitFace(); aSE.MoreFace(); aSE.NextFace()) {
					if (aSE.RejectFace(aLin) == Standard_False) {
						TopoDS_Shape aLocalShape = aSE.CurrentFace();
						TopoDS_Face CurFace = TopoDS::Face(aLocalShape);

						if (CurFace.IsSame(face)) {
							continue;
						}

						IntCurvesFace_Intersector& Intersector3d = aSE.Intersector(CurFace);
						Intersector3d.Perform(aLin, 0, RealLast());

						if (Intersector3d.IsDone()) {
							intersectNum += Intersector3d.NbPnt();
						}
					}
				}
			}
		} //end of loop on the whole solid

		if (intersectNum % 2 == 0) {
			break;
		}
	}

	bool bOut = (intersectNum / 2) % 2 == 0;
	return bOut;
}

void CuttingRailsImpl::GetOuterFaceByRayInterersect(const shape2shapeMap& map, const TopoDS_Face& benchmarkFace, const CoedgeInfo& coedgeInfos,
	TopTools_IndexedMapOfShape& visited, BRepClass3d_SolidExplorer& aSE) {

	// 这里使用 queue 实现“真正的”BFS
	std::queue<TopoDS_Face> queue;
	queue.push(benchmarkFace);

	// BFS 主循环
	while (!queue.empty())
	{
		// 取队列头
		TopoDS_Face curFace = TopoDS::Face(queue.front());
		queue.pop();

		// 若已经访问过，则跳过
		if (visited.Contains(curFace))
		{
			continue;
		}
		// 标记为已访问
		visited.Add(curFace);

		// 遍历 coedgeInfos 来找与 curFace 相邻的面
		for (auto& edge : coedgeInfos)
		{
			// 每项 edge.first 至少包含两条 coedge(面A, 面B)，edge.second 表示它们是否相反
			if (edge.first.size() != 2) {
				continue;
			}
			auto faceA = TopoDS::Face(edge.first[0]->face);
			auto faceB = TopoDS::Face(edge.first[1]->face);

			// 如果 faceA 就是当前面 curFace
			if (faceA.IsEqual(curFace) && edge.first[0]->bOuterLoop)
			{
				if (IsOuterFace(faceB, aSE)) {
					queue.push(faceB);
				}
			}

			if (faceB.IsEqual(curFace) && edge.first[1]->bOuterLoop)
			{
				if (IsOuterFace(faceA, aSE)) {
					queue.push(faceA);
				}
			}
		}
	}
}

void CuttingRailsImpl::GetFaceContinuesEdges(const TopoDS_Shape& face, const CoedgeInfo& coedgeInfos, const EdgeVectorPtr& faceCoedges) {
	(*faceCoedges).clear();
	for (auto& coedge : coedgeInfos) {
		if (coedge.first.size() == 2) {
			auto ceA = coedge.first[0];
			auto ceB = coedge.first[1];
			if (ceA->face.IsEqual(face)) {
				faceCoedges->push_back({ ceA, ceB });
			}
			else if (ceB->face.IsEqual(face)) {
				faceCoedges->push_back({ ceB, ceA });
			}
		}
		else if (coedge.first.size() == 1) {
			auto ceA = coedge.first[0];
			if (ceA->face.IsEqual(face)) {
				faceCoedges->push_back({ ceA, nullptr });
			}
		}
		else {
			for (size_t i = 0; i < coedge.first.size(); i++)
			{
				if (coedge.first[i]->face.IsEqual(face)) {
					faceCoedges->push_back({ coedge.first[i], nullptr });
				}
			}
		}
	}
}

void CuttingRailsImpl::MakeEdgeType(const CoedgeInfo& coedges) {
	int index = -1;
	for (auto& coedge : coedges) {
		index++;
		if (coedge.first.size() == 2) {
			auto ceA = coedge.first[0];
			auto ceB = coedge.first[1];
			gp_Vec vA = ceA->midNormal;
			gp_Vec vB = ceB->midNormal;
			gp_Vec crossV = vA.Crossed(vB);
			gp_Vec midV = ceA->midVec;
			if (vA.Angle(vB) < M_PI / 18) { // 10角度
				ceA->edgetType = SMOOTH;
				ceB->edgetType = SMOOTH;
				continue;
			}
			BRepBuilderAPI_MakeVertex vMaker(ceA->mid);
			auto v = vMaker.Shape();
			if (crossV.SquareMagnitude() < 1.0e-12) {
				ceA->edgetType = NOTKNOWN;
				ceB->edgetType = NOTKNOWN;
				continue;
			}

			crossV.Normalize();
			double val = midV.Dot(crossV);

			if (std::abs(val) < 0.9) {
				ceA->edgetType = NOTKNOWN;
				ceB->edgetType = NOTKNOWN;
			}
			else if (val < 0.0) {
				ceA->edgetType = AO;
				ceB->edgetType = AO;
			}
			else if (val > 0.0) {
				ceA->edgetType = TU;
				ceB->edgetType = TU;
			}
		}

		//if (coedge.first.size() > 2) {
		//	std::vector<gp_Vec> vecs;
		//	for (size_t i = 0; i < coedge.first.size(); i++)
		//	{
		//		vecs.push_back(GetNormalOfCoedgeOfFace(coedge.first[i]));
		//	}

		//	auto vA = vecs[0];
		//	bool allSmooth = true;
		//	for (size_t i = 1; i < vecs.size(); i++)
		//	{
		//		auto vB = vecs[i];
		//		if (vA.Angle(vB) > M_PI / 18) { // 10角度
		//			allSmooth = false;
		//			break;
		//		}
		//	}

		//	if (allSmooth) {
		//		for (size_t i = 0; i < coedge.first.size(); i++)
		//		{
		//			coedge.first[i]->edgetType = SMOOTH;
		//		}
		//	}

		//}
	}
}

TopoDS_Shape CuttingRailsImpl::GetVertex(const std::vector<std::pair<gp_Pnt, int>>& vertexCountA) {
	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);


	for (size_t i = 0; i < vertexCountA.size(); i++)
	{
		if (vertexCountA[i].second >= 4) {
			BRepBuilderAPI_MakeVertex vMaker(vertexCountA[i].first);
			vMaker.Build();
			BB.Add(comp, vMaker.Shape());
		}
	}

	return comp;
}

// 在shape里找相邻
void CuttingRailsImpl::GetSameClassificationPatches(const TopoDS_Shape& face, const TopTools_IndexedMapOfShape& allFaces, const Face2CoedgesMap& faceCoedgesMap,
	TopTools_IndexedMapOfShape& faceGroup, double thinDist, double tol) {


	// 这里使用 queue 实现“真正的”BFS
	std::queue<TopoDS_Shape> queue;
	queue.push(face);

	// BFS 主循环
	while (!queue.empty())
	{
		// 取队列头
		TopoDS_Face curFace = TopoDS::Face(queue.front());
		queue.pop();

		// 若已经访问过，则跳过
		if (faceGroup.Contains(curFace))
		{
			continue;
		}
		// 标记为已访问
		faceGroup.Add(curFace);

		// 遍历 coedgeInfos 来找与 curFace 相邻的面
		const auto& twinCoedges = faceCoedgesMap.FindFromKey(curFace);
		for (auto& edge : *twinCoedges)
		{
			if (!edge.second) {
				continue;
			}

			// 如果 faceA 就是当前面 curFace
			if (IsSameClassification(edge, faceCoedgesMap, thinDist, tol) && allFaces.Contains(edge.second->face)) {
				queue.push(edge.second->face);
			}
		}
	}
}

TopoDS_Shape CuttingRailsImpl::MakeSameClassificationPatchs(const TopTools_IndexedMapOfShape& shape, const Face2CoedgesMap& faceCoedgesMap,
	double thinDist, double tol, bool bConsiderLength) {

	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);

	TopTools_IndexedMapOfShape faceGroup;
	int indexStart = 1;
	for (int i = 1; i <= shape.Extent(); i++) {
		TopoDS_Face face = TopoDS::Face(shape.FindKey(i));

		if (faceGroup.Contains(face)) {
			continue;
		}

		GetSameClassificationPatches(face, shape, faceCoedgesMap, faceGroup, thinDist, tol);
		TopoDS_Shell shell;
		BB.MakeShell(shell);
		for (size_t i = indexStart; i <= faceGroup.Extent(); i++)
		{
			BB.Add(shell, faceGroup.FindKey(i));
		}
		BB.Add(comp, shell);

		indexStart = faceGroup.Extent() + 1;
	}
	return comp;
}

void CuttingRailsImpl::MakeTwoGroups(const TopoDS_Shape& patchs, TopTools_IndexedMapOfShape& outerBorderFaces, TopTools_IndexedMapOfShape& innerBorderFaces) {
	if (NbChildren(patchs) == 2) {
		TopExp_Explorer exp(patchs, TopAbs_SHELL);
		for (TopExp_Explorer expFace(exp.Current(), TopAbs_FACE); expFace.More(); expFace.Next())
		{
			TopoDS_Face curFace = TopoDS::Face(expFace.Current());
			outerBorderFaces.Add(curFace);
		}

		exp.Next();
		for (TopExp_Explorer expFace(exp.Current(), TopAbs_FACE); expFace.More(); expFace.Next())
		{
			TopoDS_Face curFace = TopoDS::Face(expFace.Current());
			innerBorderFaces.Add(curFace);
		}
	}
	else {

		BRepClass3d_SolidExplorer aSE(patchs);

		for (TopExp_Explorer exp(patchs, TopAbs_SHELL); exp.More(); exp.Next())
		{
			TopExp_Explorer expFace(exp.Current(), TopAbs_FACE);
			TopoDS_Face curFace = TopoDS::Face(expFace.Current());

			bool firstIsOuter = IsOuterFace(curFace, aSE);

			auto& targetContainer = firstIsOuter ? outerBorderFaces : innerBorderFaces;
			AddFacesToGroup(exp.Current(), targetContainer);
		}
	}
}

TopoDS_Shape CuttingRailsImpl::CollectFuzzyFaces(const TopTools_IndexedMapOfShape& fuzzyFaces,
	TopTools_IndexedMapOfShape& outerBorderFaces, TopTools_IndexedMapOfShape& innerBorderFaces,
	const Face2CoedgesMap& face2Coedges, double thinDist, double tol) {


	TopoDS_Shape fuzzyPatchs = MakeSameClassificationForFuzzyFaces(fuzzyFaces, face2Coedges, thinDist, tol);

	return fuzzyPatchs;
}

CuttingRailsImpl::AdjecentType CuttingRailsImpl::CheckFacesDistribution(const TopTools_IndexedMapOfShape& currentPatchAdjecentFaces,
	const TopTools_IndexedMapOfShape& outerBorderFaces,
	const TopTools_IndexedMapOfShape& innerBorderFaces)
{
	bool hasOuterFace = false;
	bool hasInnerFace = false;

	// 遍历 currentPatchAdjecentFaces 中的每个面
	for (Standard_Integer i = 1; i <= currentPatchAdjecentFaces.Extent(); i++)
	{
		const TopoDS_Shape& face = currentPatchAdjecentFaces(i);

		// 检查是否在外边界面中
		if (hasOuterFace || outerBorderFaces.Contains(face))
		{
			hasOuterFace = true;
		}

		// 检查是否在内边界面中
		if (hasInnerFace || innerBorderFaces.Contains(face))
		{
			hasInnerFace = true;
		}

		// 如果已经找到了同时存在于内外边界的情况，可以提前返回
		if (hasOuterFace && hasInnerFace)
		{
			return BOTH;
		}
	}

	if (hasOuterFace && !hasInnerFace) {
		return OUTER;
	}
	else if (!hasOuterFace && hasInnerFace) {
		return INNER;
	}
	else if (hasOuterFace && hasInnerFace) {
		return BOTH;
	}
	else {
		return UNKNOWNAdjecentType;
	}
}

bool CuttingRailsImpl::HasAo(const EdgeVectorPtr& edges) {
	for (const auto edge : *edges)
	{
		if (edge.first->edgetType == AO) {
			return true;
		}
	}
	return false;
}

void CuttingRailsImpl::GetMayBumpFaces(const CoedgeInfo& coedges, TopTools_IndexedMapOfShape& ret, double thinDist) {

	for (auto& coedge : coedges) {
		if (coedge.first.size() == 2) {
			auto ceA = coedge.first[0];
			auto ceB = coedge.first[1];
			if (ceA->edgetType == AO && ceA->len > 2 * thinDist) {
				auto fA = ceA->face;
				auto fB = ceB->face;
				GProp_GProps gprops;
				BRepGProp::SurfaceProperties(fA, gprops);
				double areaA = fabs(gprops.Mass());
				BRepGProp::SurfaceProperties(fB, gprops);
				double areaB = fabs(gprops.Mass());
				//面积小的是凸起
				auto f = areaA > areaB ? fB : fA;
				ret.Add(f);
			}
		}
	}
}

bool CuttingRailsImpl::IsExclusivePair(const TopoDS_Face& face1, const TopoDS_Face& face2,
	const std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair) {
	for (auto& exclusivePair : mutuallyExclusivePair) {
		if ((HasFaceInShell(face1, exclusivePair.first) && HasFaceInShell(face2, exclusivePair.second))
			|| (HasFaceInShell(face1, exclusivePair.second) && HasFaceInShell(face2, exclusivePair.first))) {
			return true;
		}
	}
	return false;
}

void CuttingRailsImpl::GetRelatePatch(const TopoDS_Face& face, const std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair,
	TopTools_IndexedMapOfShape& curPatches, TopTools_IndexedMapOfShape& otherPatches) {
	for (auto& exclusivePair : mutuallyExclusivePair) {
		if (HasFaceInShell(face, exclusivePair.first)) {
			curPatches.Add(exclusivePair.first);
			otherPatches.Add(exclusivePair.second);
		}
		else if (HasFaceInShell(face, exclusivePair.second)) {
			curPatches.Add(exclusivePair.second);
			otherPatches.Add(exclusivePair.first);
		}
	}
}

TopoDS_Shell CuttingRailsImpl::MakeFaceGroup(const TopoDS_Face& face, const Face2CoedgesMap& faceCoedgesMap,
	TopTools_IndexedMapOfShape& visitedFaceGroup,
	const TopTools_IndexedMapOfShape& innerouterFaces,
	const std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair) {
	TopoDS_Shell group;
	BRep_Builder builder;
	builder.MakeShell(group);
	// 这里使用 queue 实现“真正的”BFS
	std::queue<TopoDS_Shape> queue;
	queue.push(face);

	// BFS 主循环
	while (!queue.empty())
	{
		// 取队列头
		TopoDS_Face curFace = TopoDS::Face(queue.front());
		queue.pop();

		// 若已经访问过，则跳过
		if (visitedFaceGroup.Contains(curFace))
		{
			continue;
		}
		// 标记为已访问
		visitedFaceGroup.Add(curFace);
		builder.Add(group, curFace);
		// 遍历 coedgeInfos 来找与 curFace 相邻的面
		const auto& twinCoedges = faceCoedgesMap.FindFromKey(curFace);
		for (auto& edge : *twinCoedges)
		{
			if (!edge.second) {
				continue;
			}

			// 如果 faceA 就是当前面 curFace
			auto& face2 = TopoDS::Face(edge.second->face);
			if (innerouterFaces.Contains(face2) && !IsExclusivePair(curFace, face2, mutuallyExclusivePair))
			{
				queue.push(face2);
			}
		}
	}
	return group;
}

TopoDS_Compound CuttingRailsImpl::MakeGroups(const Face2CoedgesMap& faceCoedgesMap,
	const TopTools_IndexedMapOfShape& innerouterFaces,
	const std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair) {

	TopTools_IndexedMapOfShape visitedFaceGroup;

	TopoDS_Compound comp;
	BRep_Builder builder;
	builder.MakeCompound(comp);
	for (int i = 1; i <= innerouterFaces.Extent(); i++) {
		auto& face1 = innerouterFaces.FindKey(i);
		if (visitedFaceGroup.Contains(face1)) {
			continue;
		}
		auto& targetGroup = MakeFaceGroup(TopoDS::Face(face1), faceCoedgesMap,
			visitedFaceGroup, innerouterFaces, mutuallyExclusivePair);
		builder.Add(comp, targetGroup);
	}

	return comp;
}

void CuttingRailsImpl::MatchSeamEdges(const std::vector<ContinueEdgePtr>& unmatched, std::vector<ContinueEdgePtr>& unmatchedAfterSeamEdgeMatch, CoedgeInfo& CCMatchedMap, double tol) {

	unmatchedAfterSeamEdgeMatch.clear();

	std::set<int> usedFlag;

	for (size_t i = 0; i < unmatched.size(); i++)
	{
		if (usedFlag.find(i) != usedFlag.end()) {
			continue;
		}

		bool bMatched = false;
		const auto& current = unmatched[i];
		Handle(Geom_Surface) surf = BRep_Tool::Surface(TopoDS::Face(current->face));
		if (!surf->IsUClosed() && !surf->IsVClosed()) {
			unmatchedAfterSeamEdgeMatch.push_back(current);
			continue;
		}

		for (size_t j = i + 1; j < unmatched.size(); j++)
		{
			const auto& second = unmatched[j];
			if (!second->face.IsEqual(current->face)) {
				continue;
			}

			bool endMatch = current->start.IsEqual(second->end, tol) && current->end.IsEqual(second->start, tol);
			gp_Pnt2d secondMidUV;
			gp_Vec secondMidVec;
			bool matched = endMatch && PointDistanceToCompCurve(current->mid, second, secondMidVec, secondMidUV) < tol;

			if (matched) {
				CCMatchedMap.push_back({ { current, second }, false });
				second->UpdateMidPt(secondMidVec, secondMidUV);
				usedFlag.insert(j);
				bMatched = true;
				break;
			}
		}

		if (!bMatched) {
			unmatchedAfterSeamEdgeMatch.push_back(current);
		}
	}
}

// 辅助函数：合并face中的重复顶点
TopoDS_Face CuttingRailsImpl::MergeSameVerticesInFace(const TopoDS_Face& face) {
	double tolerance = Precision::Confusion(); // 1.0e-7

	TopoDS_Face forwardFace = face;
	forwardFace.Orientation(TopAbs_FORWARD);

	// 创建新的face
	BRep_Builder builder;
	TopoDS_Face newFace = TopoDS::Face(forwardFace.EmptyCopied());

	// 收集所有顶点
	TopTools_IndexedMapOfShape vertexMap;
	TopExp::MapShapes(forwardFace, TopAbs_VERTEX, vertexMap);

	// 创建顶点映射表
	std::map<int, int> vertexMergeMap;  // 原始顶点索引 -> 合并后顶点索引

	// 遍历顶点查找相近点
	for (int i = 1; i <= vertexMap.Extent(); i++) {
		if (vertexMergeMap.find(i) != vertexMergeMap.end()) continue;

		TopoDS_Vertex v1 = TopoDS::Vertex(vertexMap(i));
		gp_Pnt p1 = BRep_Tool::Pnt(v1);

		for (int j = i + 1; j <= vertexMap.Extent(); j++) {
			TopoDS_Vertex v2 = TopoDS::Vertex(vertexMap(j));
			gp_Pnt p2 = BRep_Tool::Pnt(v2);

			if (p1.Distance(p2) < tolerance) {
				vertexMergeMap[j] = i;
			}
		}
	}

	if (vertexMergeMap.empty()) {
		return face;
	}

	// 重建face的边界
	BRep_Builder BB;
	TopExp_Explorer wireExplorer(forwardFace, TopAbs_WIRE);
	for (; wireExplorer.More(); wireExplorer.Next()) {
		TopoDS_Wire oldWire = TopoDS::Wire(wireExplorer.Current());
		TopoDS_Wire newWire;
		builder.MakeWire(newWire);

		// 处理wire中的每条边
		TopExp_Explorer edgeExplorer(oldWire, TopAbs_EDGE);
		for (; edgeExplorer.More(); edgeExplorer.Next()) {
			TopoDS_Edge oldEdge = TopoDS::Edge(edgeExplorer.Current());
			TopoDS_Edge newEdge = oldEdge;

			// 获取边的顶点
			TopoDS_Vertex v1, v2;
			TopExp::Vertices(oldEdge, v1, v2, false);

			// 查找是否需要替换顶点
			int idx1 = vertexMap.FindIndex(v1);
			int idx2 = vertexMap.FindIndex(v2);

			if (vertexMergeMap.find(idx1) != vertexMergeMap.end() ||
				vertexMergeMap.find(idx2) != vertexMergeMap.end()) {
				// 需要创建新的边
				TopoDS_Vertex newV1 = TopoDS::Vertex(vertexMap(
					vertexMergeMap.find(idx1) != vertexMergeMap.end() ?
					vertexMergeMap[idx1] : idx1));

				TopoDS_Vertex newV2 = TopoDS::Vertex(vertexMap(
					vertexMergeMap.find(idx2) != vertexMergeMap.end() ?
					vertexMergeMap[idx2] : idx2));

				// 创建新边
				Handle(Geom_Curve) curve;
				Standard_Real f, l;
				curve = BRep_Tool::Curve(oldEdge, f, l);

				BRepBuilderAPI_MakeEdge edgeMaker(curve, newV1, newV2, f, l);
				if (edgeMaker.IsDone()) {
					newEdge = edgeMaker.Edge();
					newEdge.Orientation(oldEdge.Orientation());
				}

				Standard_Real f2d, l2d;
				Handle(Geom2d_Curve) curve2d = BRep_Tool::CurveOnSurface(oldEdge, forwardFace, f2d, l2d);
				if (!curve2d.IsNull()) {
					BB.UpdateEdge(newEdge, curve2d, newFace, Precision::Confusion());
				}
			}

			builder.Add(newWire, newEdge);
		}

		builder.Add(newFace, newWire);
	}

	newFace.Orientation(face.Orientation());

	TopTools_IndexedMapOfShape vertexMap2;
	TopExp::MapShapes(newFace, TopAbs_VERTEX, vertexMap2);

	return newFace;
}

// 原有的face处理函数
TopoDS_Face CuttingRailsImpl::MergeSameVerticesForFace(const TopoDS_Face& face) {
	Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
	if (surf->IsUClosed() || surf->IsVClosed()) {
		return MergeSameVerticesInFace(face);
	}
	else {
		return face;
	}
}

TopoDS_Compound CuttingRailsImpl::MergeSameVerticesForFaces(const TopoDS_Shape& shape) {
	// 创建新的shell
	BRep_Builder builder;
	TopoDS_Compound newShell;
	builder.MakeCompound(newShell);

	// 遍历所有face
	TopExp_Explorer explorer(shape, TopAbs_FACE);
	int i = 0;
	for (; explorer.More(); explorer.Next(), i++) {
		TopoDS_Face currentFace = TopoDS::Face(explorer.Current());

		// 处理每个face
		TopoDS_Face newFace = MergeSameVerticesForFace(currentFace);

		// 将处理后的face添加到新shell中
		builder.Add(newShell, newFace);
	}

	return newShell;
}

void CuttingRailsImpl::GetEdgesSplitInterval(const TopoDS_Shape& shape, EdgeSplitMap& edgeSplitMap)
{
	double first, last;
	for (TopExp_Explorer expFace(shape, TopAbs_FACE); expFace.More(); expFace.Next()) {
		auto& face = TopoDS::Face(expFace.Current());
		for (TopExp_Explorer expEdge(face, TopAbs_EDGE); expEdge.More(); expEdge.Next()) {
			auto& edge = TopoDS::Edge(expEdge.Current());
			Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);
			if (curve.IsNull()) {
				continue;
			}
			Handle(Geom_BSplineCurve) spline = Handle(Geom_BSplineCurve)::DownCast(curve);
			if (spline.IsNull() || spline->Continuity() != GeomAbs_C0) {
				continue;
			}

			// GeomAdaptor_Curve gac(curve);
			BRepAdaptor_Curve bac(edge);
			int nbInterval = bac.NbIntervals(GeomAbs_CN);
			if (nbInterval <= 1) {
				continue;
			}
			TColStd_Array1OfReal ts(1, nbInterval + 1);
			bac.Intervals(ts, GeomAbs_CN);
			std::vector<double> splitParams;
			splitParams.reserve(nbInterval - 1);
			for (int i = 2; i <= nbInterval; i++) {
				splitParams.push_back(ts(i));
			}
			edgeSplitMap.Add(edge, splitParams);
		}
	}
}

bool CuttingRailsImpl::bPointsEqual(const gp_Pnt& p1, const gp_Pnt& p2, const double tol) {
	// 1. 先比较x坐标
	double dx = std::abs(p1.X() - p2.X());
	if (dx > tol) {
		return false;
	}

	// 2. 再比较y坐标
	double dy = std::abs(p1.Y() - p2.Y());
	if (dy > tol) {
		return false;
	}

	// 3. 再比较z坐标
	double dz = std::abs(p1.Z() - p2.Z());
	if (dz > tol) {
		return false;
	}

	// 4. 最后检查平方距离
	double squareDistance = dx * dx + dy * dy + dz * dz;
	return squareDistance <= tol * tol;
}

void CuttingRailsImpl::GetSoftVertices(const TopoDS_Shape& shape, double tol, TopTools_IndexedMapOfShape& softVertices, EdgeSplitMap& edgeSplitMap) {
	TopTools_IndexedMapOfShape edges;
	TopExp::MapShapes(shape, TopAbs_EDGE, edges);
	double aDist;
	TopTools_IndexedDataMapOfShapeListOfShape VEMap;
	TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, VEMap);
	Standard_Real eps = tol * 0.5;
	BRepBuilderAPI_BndBoxTree aTree;
	NCollection_UBTreeFiller <Standard_Integer, Bnd_Box> aTreeFiller(aTree);
	BRepBuilderAPI_BndBoxTreeSelector aSelector;
	for (int i = 1; i <= VEMap.Extent(); i++) {
		gp_Pnt pt = BRep_Tool::Pnt(TopoDS::Vertex(VEMap.FindKey(i)));
		Bnd_Box aBox;
		aBox.Set(pt);
		aBox.Enlarge(eps);
		aTreeFiller.Add(i, aBox);
	}
	aTreeFiller.Fill();

	Handle(NCollection_BaseAllocator) allocator = NCollection_BaseAllocator::CommonBaseAllocator();
	NCollection_DataMap<TopoDS_Shape, GeomAPI_ProjectPointOnCurve*, TopTools_ShapeMapHasher> projPCMap;
	for (int i = 1; i <= edges.Extent(); i++) {

		//
		const auto& edge = TopoDS::Edge(edges(i));
		if (BRep_Tool::Degenerated(edge)) {
			continue;
		}

		double f, l;
		Handle(Geom_Curve)aC3D = BRep_Tool::Curve(edge, f, l);
		//
		const auto& pProjPC = (GeomAPI_ProjectPointOnCurve*)allocator->Allocate(sizeof(GeomAPI_ProjectPointOnCurve));
		new (pProjPC) GeomAPI_ProjectPointOnCurve();
		pProjPC->Init(aC3D, f, l);
		//
		projPCMap.Bind(edge, pProjPC);
	}

	for (int i = 1; i <= edges.Extent(); i++)
	{
		const auto& edge = TopoDS::Edge(edges(i));
		if (BRep_Tool::Degenerated(edge)) {
			continue;
		}
		Bnd_Box aGlobalBox;
		TopoDS_Vertex V1, V2;
		double first, last;
		gp_Pnt aP;

		Handle(Geom_Curve) c3d = BRep_Tool::Curve(edge, first, last);

		GeomAdaptor_Curve adptC(c3d, first, last);
		BndLib_Add3dCurve::Add(adptC, tol, aGlobalBox);
		aSelector.SetCurrent(aGlobalBox);
		aTree.Select(aSelector);
		if (!aSelector.ResInd().Extent()) continue;

		TopExp::Vertices(edge, V1, V2);
		GeomAPI_ProjectPointOnCurve* aProjector;
		projPCMap.Find(edge, aProjector);

		TColStd_ListIteratorOfListOfInteger itl(aSelector.ResInd());
		std::vector<double> splitParams;
		for (; itl.More(); itl.Next()) {
			const Standard_Integer index = itl.Value();
			const TopoDS_Vertex& vertex = TopoDS::Vertex(VEMap.FindKey(index));
			aP = BRep_Tool::Pnt(vertex);
			if (vertex.IsSame(V1) || vertex.IsSame(V2) ||
				bPointsEqual(aP, BRep_Tool::Pnt(V1), tol) || bPointsEqual(aP, BRep_Tool::Pnt(V2), tol)) {
				continue;
			}

			aProjector->Perform(aP);
			int aNbProj = aProjector->NbPoints();
			if (aNbProj) {
				aDist = aProjector->LowerDistance();
				if (aDist < tol * 1.0e-3) {
					bool bSmooth = IsSmoothAtVertex(VEMap, vertex);
					if (bSmooth) {
						softVertices.Add(vertex);
					}
					else {
						splitParams.push_back(aProjector->LowerDistanceParameter());
					}
				}
			}
		}
		if (!splitParams.empty()) {
			edgeSplitMap.Add(edge, splitParams);
		}
		aSelector.ClearResList();
	}
}

void CuttingRailsImpl::SplitEdgeByParams(const TopoDS_Edge& edge, const std::vector<double>& params, std::vector<TopoDS_Edge>& resultEdges) {
	if (edge.IsNull()) {
		return;
	}

	// 从 originalEdge 获取曲线和参数范围 [f, l]（不带 location）
	Standard_Real f = 0.0, l = 1.0;
	Handle(Geom_Curve) rawCurve = BRep_Tool::Curve(edge, f, l);
	if (rawCurve.IsNull()) {
		// 可能只存在面上的 pcurve 或无 3D 曲线。此处简单返回
		return;
	}

	Handle(Geom_Curve) curve = rawCurve;
	// 获取 Edge 的 (startV, endV)。注意到如果 orientation = Reversed，二者可能在几何上与我们直觉相反
	TopoDS_Vertex edgeStartV, edgeEndV;
	TopExp::Vertices(edge, edgeStartV, edgeEndV);

	double paramTol = Precision::PConfusion();

	// 构建最终需要的分段数组，包含 [f, l] 和用户提供的拆分点
	std::vector<Standard_Real> sortedParams = params;
	std::sort(sortedParams.begin(), sortedParams.end());
	if (sortedParams.size() > 1) {
		// 去重/过滤后的结果存放到另一个容器中
		std::vector<double> filtered;
		filtered.reserve(sortedParams.size());

		// 放入第一个参数
		filtered.push_back(sortedParams[0]);

		// 遍历剩余参数
		for (size_t i = 1; i < sortedParams.size(); ++i)
		{
			double prevVal = filtered.back();
			double currVal = sortedParams[i];
			// 如果当前值与上一个值相差小于 tol，则视为重复
			if (std::fabs(currVal - prevVal) < paramTol)
			{
				// 跳过当前值
				continue;
			}
			// 否则保留
			filtered.push_back(currVal);
		}

		// 用过滤后的覆盖原 params
		sortedParams.swap(filtered);
	}

	std::vector<Standard_Real> allParams;
	allParams.reserve(sortedParams.size() + 2);
	allParams.push_back(f);
	
	for (auto p : sortedParams) {
		if (p > f + paramTol && p < l - paramTol) {
			allParams.push_back(p);
		}
	}
	allParams.push_back(l);

	// 构造子 Edge
	BRep_Builder builder;

	// 当前段的起点顶点
	TopoDS_Vertex currStartV = edgeStartV;
	bool sameParameter = BRep_Tool::SameParameter(edge);
	bool sameRange = BRep_Tool::SameRange(edge);
	double edgeTol = BRep_Tool::Tolerance(edge);
	bool edgeReversed = edge.Orientation() == TopAbs_REVERSED;
	for (size_t i = 0; i < allParams.size() - 1; ++i)
	{
		Standard_Real segStart = allParams[i];
		Standard_Real segEnd = allParams[i + 1];

		// 构造裁剪曲线
		Handle(Geom_TrimmedCurve) subCurve =
			new Geom_TrimmedCurve(curve, segStart, segEnd);

		// 构造子 Edge
		TopoDS_Edge subEdge = TopoDS::Edge(edge.EmptyCopied());

		builder.Range(subEdge, segStart, segEnd);

		// 计算终点，即 segEnd 处的 3D 
		gp_Pnt startPt = curve->Value(segStart);
		gp_Pnt endPt = curve->Value(segEnd);
		// 如果是最后一段，把终点设为原始 Edge 的终点；否则创建一个新的顶点
		TopoDS_Vertex currEndV;
		if (i == allParams.size() - 2) {
			currEndV = edgeEndV;
		}
		else {
			// 新建顶点（不带 location，后面 Edge 自己有 location）
			builder.MakeVertex(currEndV, endPt, Precision::Confusion());
		}
		// 将起点顶点、终点顶点加到子 Edge
		if (edgeReversed) {
			builder.Add(subEdge, currStartV.Oriented(TopAbs_REVERSED));
			builder.Add(subEdge, currEndV.Oriented(TopAbs_FORWARD));
		}
		else {
			builder.Add(subEdge, currStartV.Oriented(TopAbs_FORWARD));
			builder.Add(subEdge, currEndV.Oriented(TopAbs_REVERSED));
		}
		

		// 检查
		//double first, last;
		//Handle(Geom_Curve) subC = BRep_Tool::Curve(subEdge, first, last);
		//BRepMesh_IncrementalMesh mesher(subEdge, 0.05, Standard_False, 0.05);
		//// 获取多边形化后的数据
		//TopLoc_Location loc;
		//Handle(Poly_Polygon3D) polygon = BRep_Tool::Polygon3D(subEdge, loc);
		//std::vector<gp_Pnt> discretePnts;
		//const TColgp_Array1OfPnt& nodes = polygon->Nodes();
		//for (int i = 1; i <= polygon->NbNodes(); i++) {
		//	discretePnts.push_back(nodes(i));
		//}

		// 收集结果
		resultEdges.push_back(subEdge);
		

		// 下一段的起点顶点
		currStartV = currEndV;
	}

	if (edgeReversed) {
		std::reverse(resultEdges.begin(), resultEdges.end());
	}
	return;
}

TopoDS_Face CuttingRailsImpl::SplitFaceEdges(const TopoDS_Face& face, const EdgeSplitMap& edgeSplitParams) {
	TopTools_IndexedMapOfShape faceEdges;
	TopExp::MapShapes(face, TopAbs_EDGE, faceEdges);
	bool needSplit = false;
	for (int i = 1; i <= edgeSplitParams.Extent(); i++) {
		if (faceEdges.Contains(edgeSplitParams.FindKey(i))) {
			needSplit = true;
			break;
		}
	}

	if (!needSplit) {
		return face;
	}

	BRep_Builder builder;
	// 用于最终输出的新面
	TopoDS_Face newFace = TopoDS::Face(face.EmptyCopied());
	bool isFaceReversed = newFace.Orientation() == TopAbs_REVERSED;
	// 遍历原始面上的所有Edge
	for (TopoDS_Iterator expWire(face, false, false); expWire.More(); expWire.Next())
	{
		if (expWire.Value().ShapeType() != TopAbs_WIRE) {
			continue;
		}
		auto& curWire = expWire.Value();
		TopoDS_Shape newWire = curWire.EmptyCopied();
		bool isWireReversed = newWire.Orientation() == TopAbs_REVERSED;
		for (TopoDS_Iterator expEdge(curWire, false, false); expEdge.More(); expEdge.Next()) {
			if (expEdge.Value().ShapeType() != TopAbs_EDGE) {
				continue;
			}
			TopoDS_Edge originalEdge = TopoDS::Edge(expEdge.Value());
			// 检查该Edge是否在映射中
			if (!edgeSplitParams.Contains(originalEdge))
			{
				// 如果不在映射中，则无需拆分，直接加入 newWire
				builder.Add(newWire, isWireReversed ? originalEdge.Reversed() : originalEdge);
				continue;
			}

			// 如果 Edge 在映射中，则根据打断点进行拆分
			const std::vector<double>& params = edgeSplitParams.FindFromKey(originalEdge);
			std::vector<TopoDS_Edge> subEdges;
			SplitEdgeByParams(originalEdge, params, subEdges);
			for (auto& subEdge : subEdges) {
				builder.Add(newWire, isWireReversed ? subEdge.Reversed() : subEdge);
			}
		}
		// 将完整的Wire附加到新的Face
		builder.Add(newFace, isFaceReversed ? newWire.Reversed() : newWire);
	}

	return newFace;
}

TopoDS_Shape CuttingRailsImpl::SplitShapeEdges(const TopoDS_Shape& shape, const EdgeSplitMap& splitParams) {
	if (splitParams.Size() == 0) {
		return shape;
	}

	TopoDS_Shape newShape = shape.EmptyCopied();
	BRep_Builder builder;
	TopoDS_Iterator iterator(shape, false, false);
	int i = 1;
	bool isShapeReversed = shape.Orientation() == TopAbs_REVERSED;
	for (; iterator.More(); iterator.Next(), i++) {
		auto& curShape = iterator.Value();
		if (curShape.ShapeType() != TopAbs_FACE) {
			continue;
		}

		auto& newFace = SplitFaceEdges(TopoDS::Face(curShape), splitParams);
		builder.Add(newShape, isShapeReversed ? newFace.Reversed() : newFace);
	}
	return newShape;
}

CuttingRailsImpl::ContinueEdgePtr CuttingRailsImpl::GetTwinContinueEdge(const ContinueEdgePtr& continueEdgePtr, const CoedgeInfo& edges) {
	for (const auto& edge : edges)
	{
		if (edge.first.size() == 2) {
			auto ceA = edge.first[0];
			auto ceB = edge.first[1];
			if (ceA == continueEdgePtr) {
				return ceB;
			}

			if (ceB == continueEdgePtr) {
				return ceA;
			}
		}
	}

	return nullptr;
}

TopoDS_Shape CuttingRailsImpl::GetBorderFacesForExtrudeShape(const TopoDS_Shape& shape) {
	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);

	gp_Vec vec;
	bool bFound = GetExtrudeDirection(shape, vec);
	if (bFound) {
		gp_Dir extrudeDir(vec);  // 将挤出方向转换为单位向量

		for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
			gp_Dir normalStart, normalMiddle, normalEnd;
			ComputeFaceNormals(exp.Current(), normalStart, normalMiddle, normalEnd);

			// 计算每个法向与挤出方向的点积的绝对值
			// 如果接近0，说明它们垂直
			Standard_Real dotStart = Abs(normalStart.Dot(extrudeDir));
			Standard_Real dotMiddle = Abs(normalMiddle.Dot(extrudeDir));
			Standard_Real dotEnd = Abs(normalEnd.Dot(extrudeDir));

			// 设定一个容差值，判断是否垂直
			const Standard_Real angleTolerance = 0.1;  // 6度

			// 只有当所有点的法向都垂直于挤出方向时，才添加这个面
			if ((dotStart <= angleTolerance &&
				dotMiddle <= angleTolerance &&
				dotEnd <= angleTolerance)) {
				BB.Add(comp, exp.Current());
			}
		}
	}

	return comp;
}

TopoDS_Compound CuttingRailsImpl::showEdgesType(const CoedgeInfo& coedgeInfoVector) {

	BRep_Builder BB;

	TopoDS_Compound comp;
	TopoDS_Compound unmatched;
	TopoDS_Compound ao;
	TopoDS_Compound tu;
	TopoDS_Compound smooth;
	TopoDS_Compound unknown;

	BB.MakeCompound(comp);
	BB.MakeCompound(unmatched);
	BB.MakeCompound(ao);
	BB.MakeCompound(tu);
	BB.MakeCompound(smooth);
	BB.MakeCompound(unknown);


	for (size_t i = 0; i < coedgeInfoVector.size(); i++)
	{
		if (coedgeInfoVector[i].first.size() == 1) {
			BB.Add(unknown, coedgeInfoVector[i].first[0]->compedges);
		}

		if (coedgeInfoVector[i].first.size() == 2) {
			auto tp = coedgeInfoVector[i].first[0]->edgetType;
			if (tp == AO) {
				BB.Add(ao, coedgeInfoVector[i].first[0]->compedges);
			}
			else if (tp == TU) {
				BB.Add(tu, coedgeInfoVector[i].first[0]->compedges);
			}
			else if (tp == SMOOTH) {
				BB.Add(smooth, coedgeInfoVector[i].first[0]->compedges);
			}
			else {
				BB.Add(unknown, coedgeInfoVector[i].first[0]->compedges);
			}
		}
	}

	BB.Add(comp, unmatched);
	BB.Add(comp, ao);
	BB.Add(comp, tu);
	BB.Add(comp, smooth);
	BB.Add(comp, unknown);

	return comp;
}

void CuttingRailsImpl::FixEdgeInfo(const CoedgeInfo& coedgeInfoVector, const TopTools_IndexedMapOfShape& reversedFaces,
	const TopTools_IndexedDataMapOfShapeShape& faceMapping,
	CoedgeInfo& coedgeInfoVectorAfterFaceOritationAdjust) {
	coedgeInfoVectorAfterFaceOritationAdjust.clear();
	coedgeInfoVectorAfterFaceOritationAdjust.resize(coedgeInfoVector.size());

	//更新反向的面和coedges里的面
	TopTools_IndexedMapOfShape reversedFacesAfterMapping;
	for (size_t i = 1; i <= reversedFaces.Extent(); i++)
	{
		reversedFacesAfterMapping.Add(faceMapping.FindFromKey(reversedFaces.FindKey(i)));
	}

	for (size_t i = 0; i < coedgeInfoVector.size(); i++)
	{
		const std::vector<ContinueEdgePtr>& coedges = coedgeInfoVector[i].first;
		for (const auto& coedge : coedges) {
			if (faceMapping.Contains(coedge->face)) {
				coedge->face = faceMapping.FindFromKey(coedge->face);
			}
		}
	}
	//完成更新

	for (size_t i = 0; i < coedgeInfoVector.size(); i++)
	{
		const std::vector<ContinueEdgePtr>& coedges = coedgeInfoVector[i].first;
		if (coedges.size() == 2) {

			bool reverse1 = reversedFacesAfterMapping.Contains(coedges[0]->face);
			bool reverse2 = reversedFacesAfterMapping.Contains(coedges[1]->face);
			if (reverse1 && reverse2) {
				coedgeInfoVectorAfterFaceOritationAdjust[i] = { { coedges[0]->Reversed(), coedges[1]->Reversed() }, coedgeInfoVector[i].second };
			}
			else if (reverse1) {
				coedgeInfoVectorAfterFaceOritationAdjust[i] = { { coedges[0]->Reversed(), coedges[1] }, !coedgeInfoVector[i].second };
			}
			else if (reverse2) {
				coedgeInfoVectorAfterFaceOritationAdjust[i] = { { coedges[0], coedges[1]->Reversed() }, !coedgeInfoVector[i].second };
			}
			else {
				coedgeInfoVectorAfterFaceOritationAdjust[i] = coedgeInfoVector[i];
			}
			coedgeInfoVectorAfterFaceOritationAdjust[i].first[0]->twinface = coedgeInfoVectorAfterFaceOritationAdjust[i].first[1]->face;
			coedgeInfoVectorAfterFaceOritationAdjust[i].first[1]->twinface = coedgeInfoVectorAfterFaceOritationAdjust[i].first[0]->face;
		}
		else {
			// 更新reverse
			std::vector<ContinueEdgePtr> newCoedges;
			for (const auto& coedge : coedges) {
				if (reversedFacesAfterMapping.Contains(coedge->face)) {
					newCoedges.push_back(coedge->Reversed());
				}
				else {
					newCoedges.push_back(coedge);
				}
			}
			coedgeInfoVectorAfterFaceOritationAdjust[i] = { newCoedges, coedgeInfoVector[i].second };
		}
	}
}

void CuttingRailsImpl::GetFace2ContinueEdges(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfoVector, Face2CoedgesMap& face2Coedges) {
	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
		TopoDS_Shape face = exp.Current();
		EdgeVectorPtr coedges = std::make_shared<EdgeVector>();
		GetFaceContinuesEdges(face, coedgeInfoVector, coedges);
		face2Coedges.Add(face, coedges);
	}
}

void CuttingRailsImpl::GetFace2ContinueEdges(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfoVector, Face2CoedgesMap& face2Coedges, TopTools_IndexedMapOfShape& allFaces) {
	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
		TopoDS_Shape face = exp.Current();
		EdgeVectorPtr coedges = std::make_shared<EdgeVector>();
		GetFaceContinuesEdges(face, coedgeInfoVector, coedges);
		face2Coedges.Add(face, coedges);
		allFaces.Add(face);
	}
}

TopTools_IndexedMapOfShape CuttingRailsImpl::GetReversedFaces(const TopoDS_Shape& shape, const TopTools_IndexedMapOfShape& faces) {
	TopTools_IndexedMapOfShape reversedFaces;
	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
		if (!faces.Contains(exp.Current())) {
			reversedFaces.Add(exp.Current());
		}
	}
	return reversedFaces;
}

TopoDS_Shape CuttingRailsImpl::showAllContinuesEdges(const std::vector<ContinueEdgePtr>& allWireContinuesEdges) {
	BRep_Builder BB;
	TopoDS_Compound comp;
	BB.MakeCompound(comp);

	for (size_t i = 0; i < allWireContinuesEdges.size(); i++)
	{
		BB.Add(comp, allWireContinuesEdges[i]->compedges);
	}
	return comp;
}

void CuttingRailsImpl::GetCuttingRails(const TopoDS_Shape& shape, TopoDS_Shape& outerBoundFaces, TopoDS_Shape& innerBoundFaces, TopoDS_Compound& unMatchedEdges,
	std::vector<std::vector<ContinueEdgePtr>>& loops, std::vector<std::vector<ContinueEdgePtr>>& twinLoops, const CuttingRailsOption& option, double tol) {
	BRep_Builder BB;
	TopExp_Explorer explorer;
	TopTools_IndexedMapOfShape cuttingFaces;
	TopTools_IndexedMapOfShape borderFaces;
	TopTools_IndexedMapOfShape outerBorderFaces;
	TopTools_IndexedMapOfShape innerBorderFaces;
	TopTools_IndexedMapOfShape fuzzyFaces;
	TopTools_IndexedMapOfShape bumpFaces;
	TopTools_IndexedMapOfShape allFaces;
	//TopTools_IndexedDataMapOfShapeListOfShape EFMap;
	//TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, EFMap);
	CoedgeInfo coedgeInfoVector;
	std::vector<double> sortedFaceAreas;
	TopoDS_Shape fixedShape = SortFaceByArea(shape, sortedFaceAreas); // 获取最大的面
	// BRepTools_Write("E:/shape.brep", (void*)&fixedShape);
	fixedShape = MergeSameVerticesForFaces(fixedShape);

	//begin data repair
	TopTools_IndexedDataMapOfShapeShape faceMapping;
	fixedShape = FixFaceWiresOrder(fixedShape, faceMapping);
	std::vector<ContinueEdgePtr> allWireContinuesEdges;

	// 对c0连续的nurbs 做打断
	/*EdgeSplitMap splineSplitParams;
	GetEdgesSplitInterval(fixedShape, splineSplitParams);
	fixedShape = SplitFaces(fixedShape, splineSplitParams);*/

	std::vector<ContinueEdgePtr> firstUnmatched;
	std::vector<ContinueEdgePtr> unmatchedAfterSeamEdge;
	std::vector<ContinueEdgePtr> unmatchedAfterMergeEdge;

	TopTools_IndexedMapOfShape softVertices;
	EdgeSplitMap edgeSplitParams;
	GetSoftVertices(fixedShape, tol, softVertices, edgeSplitParams);
	fixedShape = SplitShapeEdges(fixedShape, edgeSplitParams);
	//TopoDS_Shape sofeVertices = showTopos(softVertices);
	std::shared_ptr<KDTree> tree = MakeKDTree(fixedShape);
	MakeContinueEdges(tree, fixedShape, allWireContinuesEdges, softVertices, tol);
	//auto allEdges = showAllContinuesEdges(allWireContinuesEdges);
	MatchCoedges(tree, allWireContinuesEdges, coedgeInfoVector, firstUnmatched, tol, unMatchedEdges);
	if (!firstUnmatched.empty()) {
		MatchSeamEdges(firstUnmatched, unmatchedAfterSeamEdge, coedgeInfoVector, tol);
		// MatchMergedEdges(tree, unmatchedAfterSeamEdge, coedgeInfoVector, unmatchedAfterMergeEdge, tol, unMatchedEdges);
		for (auto& unmatched : unmatchedAfterSeamEdge) {
			coedgeInfoVector.push_back({ { unmatched }, false });
		}
	}

	double thinDist = ComputeThicknessByRay(fixedShape, TopoDS::Face(GetFirstSubShape(fixedShape)), coedgeInfoVector);
	// std::cout << "thin: " << thinDist << std::endl;
	if (thinDist < 0) {
		thinDist = 0.0;
	}
	else {
		if (thinDist > 30 || option.singleBorderType != SingleBorderType::None) {
			thinDist = 0;
		}
		else {
			tol = std::min(tol, 0.95 * thinDist);
		}
	}

	TopTools_IndexedMapOfShape reversedFaces;
	fixedShape = FixFaceOritation(fixedShape, coedgeInfoVector, reversedFaces);

	gp_Vec extrudeDir;
	bool hasExtrudeDir = GetSortedShapeExtrudeDir(fixedShape, sortedFaceAreas, extrudeDir);
	TopoDS_Solid solid;
	BB.MakeSolid(solid);
	BB.Add(solid, fixedShape);
	bool isSingle = thinDist == 0 && hasExtrudeDir;
	int state = SolidInOutStatus(solid, isSingle, extrudeDir);
	bool reverse = false;
	if (isSingle) {
		// 单壁处理，考虑singleBorderType
		// 当单壁类型未知时，当做外壁处理
		// 当单壁类型为外壁且朝向为内时，翻转
		// 当单壁类型为内壁且朝向为外时，翻转
		reverse = (state == TopAbs_IN && option.singleBorderType != SingleBorderType::Inner) || (state == TopAbs_OUT && option.singleBorderType == SingleBorderType::Inner);
	}
	else {
		// 面向封闭实体模型处理
		reverse = state == TopAbs_IN;
	}
	if (reverse) {
		fixedShape.Reverse();
		reversedFaces = GetReversedFaces(fixedShape, reversedFaces);
	}

	fixedShape = FixFaceWiresOrder(fixedShape, faceMapping);
	// end data repair
	// BRepTools_Write("E:/shape.brep", (void*)&fixedShape);
	CoedgeInfo coedgeInfoVectorAfterFaceOritationAdjust;
	FixEdgeInfo(coedgeInfoVector, reversedFaces, faceMapping, coedgeInfoVectorAfterFaceOritationAdjust); //修正方向, 同时补充twin face

	Face2CoedgesMap face2Coedges;
	GetFace2ContinueEdges(fixedShape, coedgeInfoVectorAfterFaceOritationAdjust, face2Coedges, allFaces);

	MakeEdgeType(coedgeInfoVectorAfterFaceOritationAdjust);

	//TopoDS_Compound types = showEdgesType(coedgeInfoVectorAfterFaceOritationAdjust);// unmatch ao tu smooth unknown

	fixedShape = MakeSameClassificationPatchs(allFaces, face2Coedges, thinDist, tol, true);

	fixedShape = SortShellByArea(fixedShape);

	std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>> mutuallyExclusivePair;

	GetInnerOuterShellByRatio(fixedShape, face2Coedges, borderFaces,
		cuttingFaces, fuzzyFaces, mutuallyExclusivePair, option.bExtrudeShape, thinDist, tol);

	if (thinDist == 0.0) {
		mutuallyExclusivePair.clear();
	}

	TopoDS_Compound comp = MakeGroups(face2Coedges, borderFaces,
		mutuallyExclusivePair);
	ClassifyInOutGroups(comp, outerBoundFaces, innerBoundFaces, hasExtrudeDir, extrudeDir, option.singleBorderType);

	std::vector<ContinueEdgePtr> candidateRails = GetCandidateRails(
		option.singleBorderType == SingleBorderType::Inner ? innerBoundFaces : outerBoundFaces, face2Coedges, option.bNeedCuttingFaces);

	std::vector<ContinueEdgePtr> remains = FindLoops(candidateRails, loops, tol);

	twinLoops.clear();
	twinLoops.reserve(loops.size());
	for (const auto& loop : loops)
	{
		std::vector<ContinueEdgePtr> twinLoop;
		twinLoop.reserve(loop.size());
		for (auto it = loop.rbegin(); it != loop.rend(); ++it) {
			const auto& continueEdge = *it;
			const auto& twin = GetTwinContinueEdge(continueEdge, coedgeInfoVectorAfterFaceOritationAdjust);
			twinLoop.push_back(twin);
		}

		twinLoops.push_back(twinLoop);
	}

	return;
}

void CuttingRailsImpl::DiscretizeCurve(const Handle(Geom_Curve)& curve, double start, double end, int numPoints, std::vector<PointWithParam>& rets) {
	std::vector<double> parameterValues;

	int segments = numPoints - 1;
	// 添加起始点的参数
	parameterValues.push_back(start);

	// 遍历计算参数值
	for (int i = 1; i <= numPoints - 1; i++) {
		double t = start + (static_cast<double>(i) / segments) * (end - start);

		// 获取当前点的切向量
		gp_Vec currentTangent;
		curve->D1(t, gp_Pnt(), currentTangent);
		currentTangent.Normalize();

		// 获取前一个点的切向量
		gp_Vec prevTangent;
		curve->D1(parameterValues.back(), gp_Pnt(), prevTangent);
		prevTangent.Normalize();

		// 计算角度
		double dot = currentTangent.Dot(prevTangent);

		// 如果角度大于阈值，添加该点
		if (dot <= 0.995) {// cos(0.1) ≈ 0.995
			parameterValues.push_back(t);
		}
	}

	// 添加终点的参数
	parameterValues.push_back(end);

	// 计算所有参数值对应的点坐标
	for (double t : parameterValues) {
		rets.push_back(PointWithParam(t, curve->Value(t)));
	}

	return;
}

void CuttingRailsImpl::GetEdgePointsAndNormals(const TopoDS_Edge& edge,
	const TopoDS_Face& face,
	const int discreteNum,
	std::vector<gp_Pnt>& pnts, std::vector<gp_Dir>& vecs)
{
	// 1. 获取edge的几何信息
	Standard_Real first, last;
	Handle(Geom2d_Curve) pcurve = BRep_Tool::CurveOnSurface(edge, face, first, last);
	Handle(Geom2d_BoundedCurve) trimmed = new Geom2d_TrimmedCurve(pcurve, first, last);
	Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
	Handle(Geom_CurveOnSurface) curveOnsurface = new Geom_CurveOnSurface(trimmed, surface);

	//离散
	//const Handle(Geom_Curve)& curve, double start, double end, int numPoints, std::vector<PointWithParam>& rets
	std::vector<PointWithParam> pointsWithParams;
	DiscretizeCurve(curveOnsurface, first, last, discreteNum, pointsWithParams);

	for (const auto& p : pointsWithParams)
	{
		pnts.push_back(p.point);
		gp_Pnt2d pt2d;
		pcurve->D0(p.t, pt2d);
		gp_Dir aDN;
		FaceNormal(face, pt2d.X(), pt2d.Y(), aDN);
		vecs.push_back(aDN);
	}
	return;
}

void CuttingRailsImpl::GetPointsAndNormalFromWires(std::vector<std::vector<ContinueEdgePtr>>& wires,
	const int discreteNum,
	std::vector<std::vector<gp_Pnt>>& pnts, std::vector<std::vector<gp_Dir>>& normals) {

	for (const auto& wire : wires) {
		std::vector<gp_Pnt> wirePnts;
		std::vector<gp_Dir> wireNormals;
		for (const auto& continueEdge : wire) {
			for (const auto edge : continueEdge->edges)
			{
				GetEdgePointsAndNormals(edge, TopoDS::Face(continueEdge->face), discreteNum, wirePnts, wireNormals);
			}
		}
		pnts.push_back(wirePnts);
		normals.push_back(wireNormals);
	}
}

bool CuttingRailsImpl::HasFace(const TopoDS_Shape& shape) {
	if (shape.IsNull()) return false;

	TopExp_Explorer exp(shape, TopAbs_FACE);
	return exp.More();  // 如果有face，More()返回true
}

bool CuttingRailsImpl::FindCommonDirection(const TopoDS_Shape& shape, gp_Vec& dir,
	const Standard_Real theAngularTolerance) {
	std::vector<gp_Vec> candidates;
	bool firstFace = true;
	CuttingRailsImpl util;

	for (TopExp_Explorer faceExp(shape, TopAbs_FACE); faceExp.More(); faceExp.Next()) {
		TopoDS_Face face = TopoDS::Face(faceExp.Current());
		gp_Vec dir1, dir2;

		// 获取面的挤压方向
		if (!util.IsExtrudeFace(face, dir1, dir2, theAngularTolerance)) {
			return false;
		}

		if (firstFace) {
			// 处理第一个面
			candidates.push_back(dir1);

			// 如果两个方向不同，都加入候选
			if (!dir1.XYZ().IsEqual(dir2.XYZ(), 1.0e-12)) {
				candidates.push_back(dir2);
			}
			firstFace = false;
		}
		else {
			// 处理后续面
			std::vector<gp_Vec> newCandidates;

			// 筛选与当前面方向匹配的候选
			for (const auto& candidate : candidates) {
				if (candidate.IsParallel(dir1, theAngularTolerance) ||
					candidate.IsParallel(dir2, theAngularTolerance)) {
					newCandidates.push_back(candidate);
				}
			}

			candidates = std::move(newCandidates);

			// 如果没有共同方向，返回false
			if (candidates.empty()) {
				return false;
			}

			// 如果只剩一个候选方向，提前结束
			if (candidates.size() == 1) {
				break;
			}
		}
	}

	// 如果有候选方向，取第一个作为结果
	if (!candidates.empty()) {
		dir.SetXYZ(candidates.front().XYZ());
		return true;
	}

	return false;
}

bool CuttingRailsImpl::FindDirectionBySurfaceNormalCross(const TopoDS_Shape& shape, gp_Vec& direction)
{
	TopExp_Explorer faceExp(shape, TopAbs_FACE);
	if (!faceExp.More()) {
		return false;
	}

	// 获取第一个面的法向
	TopoDS_Face firstFace = TopoDS::Face(faceExp.Current());
	Handle(Geom_Surface) firstSurface = BRep_Tool::Surface(firstFace);
	Standard_Real uMin1, uMax1, vMin1, vMax1;
	BRepTools::UVBounds(firstFace, uMin1, uMax1, vMin1, vMax1);

	// 计算第一个面中点的法向
	Standard_Real uMid1 = (uMin1 + uMax1) * 0.5;
	Standard_Real vMid1 = (vMin1 + vMax1) * 0.5;
	gp_Pnt point1;
	gp_Vec du1, dv1;
	firstSurface->D1(uMid1, vMid1, point1, du1, dv1);
	gp_Vec normal1 = du1.Crossed(dv1);
	if (normal1.Magnitude() < Precision::Confusion()) {
		return false;
	}
	normal1.Normalize();

	// 前进到第二个面
	faceExp.Next();

	// dot 阈值（对应夹角为30°）
	const double dotThreshold = std::cos(M_PI / 6.0); // ≈ 0.8660254

	// 遍历其余的面
	for (; faceExp.More(); faceExp.Next()) {
		TopoDS_Face face = TopoDS::Face(faceExp.Current());
		Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
		Standard_Real uMin, uMax, vMin, vMax;
		BRepTools::UVBounds(face, uMin, uMax, vMin, vMax);

		// 计算当前面中点的法向
		Standard_Real uMid = (uMin + uMax) * 0.5;
		Standard_Real vMid = (vMin + vMax) * 0.5;
		gp_Pnt point;
		gp_Vec du, dv;
		surface->D1(uMid, vMid, point, du, dv);
		gp_Vec normal2 = du.Crossed(dv);
		if (normal2.Magnitude() < Precision::Confusion()) {
			continue;
		}
		normal2.Normalize();

		// 计算两个法向的点乘
		double dotValue = normal1.Dot(normal2);

		// 如果点乘小于阈值，说明夹角大于30度，立即计算 resultDir 并返回
		if (dotValue < dotThreshold) {
			gp_Vec resultDir = normal1.Crossed(normal2);
			if (resultDir.Magnitude() > Precision::Confusion()) {
				resultDir.Normalize();
				direction.SetXYZ(resultDir.XYZ());
				return true;
			}
		}
	}

	return false;
}

//要求面积排过序
bool CuttingRailsImpl::GetExtrudeDirection(const TopoDS_Shape& shape, gp_Vec& dir) {

	double theAngularTolerance = 0.05;
	bool suc = FindCommonDirection(shape, dir, theAngularTolerance);
	/*if (!suc) {
		suc = FindDirectionBySurfaceNormalCross(shape, dir);
	}*/
	return suc;
}

TopoDS_Compound CuttingRailsImpl::copyShape(const TopoDS_Shape& shape) {
	// 创建一个新的 Compound，用来容纳复制后的面
	BRep_Builder builder;
	TopoDS_Compound comp;
	builder.MakeCompound(comp);

	// 遍历原 Shell 的所有面
	TopTools_MapOfShape faces;
	for (TopExp_Explorer faceExp(shape, TopAbs_FACE); faceExp.More(); faceExp.Next())
	{
		// 获取当前的面
		TopoDS_Face face = TopoDS::Face(faceExp.Current());

		if (faces.Contains(face)) {
			continue;
		}
		// 使用 BRepBuilderAPI_Copy 对 face 进行复制
		BRepBuilderAPI_Copy copyTool(face);
		TopoDS_Shape copiedFace = copyTool.Shape();

		// 将复制得到的独立面加入到 Compound 中
		builder.Add(comp, copiedFace);
		faces.Add(face);
	}

	return comp;
}

bool CuttingRailsImpl::IsCircularCurve(const std::vector<gp_Pnt>& points, gp_Pnt& center, gp_Vec& normal, double squareTol) {
	// 取三个点确定一个平面和可能的圆心
	int numPoints = static_cast<int>(points.size());
	gp_Pnt p1(points[0]);
	gp_Pnt p2(points[numPoints / 4]);
	gp_Pnt p3(points[numPoints / 2]);
	double radius = 0;

	GC_MakeCircle mc(p1, p2, p3);
	if (mc.IsDone()) {
		Handle(Geom_Circle) circle = mc.Value();
		center = circle->Location();
		normal = circle->Axis().Direction();
		radius = circle->Radius();
	}
	else {
		return false;
	}
	gp_Pln plane(center, normal);
	double squareRaidus = radius * radius;
	// 验证所有点是否在圆上
	for (const auto& p : points) {
		double dist = p.SquareDistance(center);
		if (std::abs(dist - squareRaidus) > squareTol)
			return false;

		// 检查点是否在平面上
		if (std::abs(plane.SquareDistance(p)) > squareTol)
			return false;
	}

	return true;
}

int CuttingRailsImpl::IsPlanerCurve(const std::vector<gp_Pnt>& points, gp_Vec& normal, double squareTol)
{
	if (points.size() < 3) {
		return 0; // 少于3个点无法确定平面
	}

	// 寻找前三个不共线的点来确定平面
	gp_Pnt p1 = points[0];
	gp_Pnt p2;
	gp_Pnt p3;

	// 找第二个点（与第一个点不重合）
	int idx2 = -1;
	for (int i = 1; i < points.size(); ++i) {
		if (p1.SquareDistance(points[i]) > squareTol) {
			p2 = points[i];
			idx2 = i;
			break;
		}
	}
	if (idx2 == -1) {
		return 0; // 所有点都重合
	}

	// 找第三个点（与前两个点不共线）
	int idx3 = -1;
	gp_Vec v12(p1, p2);
	v12.Normalize();
	for (int i = idx2 + 1; i < points.size(); ++i) {
		gp_Vec v13(p1, points[i]);
		if (v13.SquareMagnitude() <= gp::Resolution()) {
			continue;
		}
		v13.Normalize();
		gp_Vec cross = v12.Crossed(v13);
		if (cross.SquareMagnitude() > 1.0e-12) {
			p3 = points[i];
			idx3 = i;
			normal = cross.Normalized();
			break;
		}
	}
	if (idx3 == -1) {
		return 0; // 所有点共线
	}

	// 检查其余所有点是否都在这个平面上
	for (int i = 0; i < points.size(); ++i) {
		if (i == 0 || i == idx2 || i == idx3) {
			continue; // 跳过用于确定平面的三个点
		}

		// 计算点到平面的距离的平方
		gp_Vec v1p(p1, points[i]);
		if (v1p.SquareMagnitude() <= gp::Resolution()) {
			continue;
		}
		v1p.Normalize();
		double dot = v1p.Dot(normal);
		double distSquare = dot * dot;

		if (distSquare > squareTol) {
			return -1; // 点不在平面上
		}
	}

	return 1;
}

bool CuttingRailsImpl::IsCircularCurve(const Handle(Geom_Curve)& curve, double uMin, double uMax, gp_Pnt& center,
	double squareTol) {
	// 取曲线上的多个点进行判断
	int numPoints = 16;
	std::vector<gp_Pnt> points;
	double du = (uMax - uMin) / (numPoints - 1);

	// 采样点
	for (int i = 0; i < numPoints; i++) {
		double u = uMin + i * du;
		gp_Pnt p = curve->Value(u);
		points.push_back(p);
	}

	gp_Vec normal;
	return IsCircularCurve(points, center, normal, squareTol);
}

int CuttingRailsImpl::IsPlanerCurve(const Handle(Geom_Curve)& curve, double uMin, double uMax, gp_Vec& normal,
	double squareTol) {
	// 取曲线上的多个点进行判断
	int numPoints = 16;
	std::vector<gp_Pnt> points;
	double du = (uMax - uMin) / (numPoints - 1);

	// 采样点
	for (int i = 0; i < numPoints; i++) {
		double u = uMin + i * du;
		gp_Pnt p = curve->Value(u);
		points.push_back(p);
	}

	return IsPlanerCurve(points, normal, squareTol);
}

bool CuttingRailsImpl::SampleContinue(const ContinueEdgePtr& continueEdge, int sampleNum, std::vector<gp_Pnt>& samplePnts) {
	samplePnts.clear();
	samplePnts.reserve(sampleNum * continueEdge->edges.size());
	double delta = continueEdge->len / (sampleNum - 1);
	double len = 0;

	for (int i = 0, iLen = continueEdge->edges.size(); i < iLen; i++) {
		double curLen = continueEdge->increaseLen[i];
		if (curLen < len) {
			continue; // 跳过
		}
		// 计算当前edge的采样数量
		int sampleCount = 0;
		while (len < curLen) {
			len += delta;
			sampleCount++;
		}
		TopoDS_Edge edge = continueEdge->edges[i];
		// 计算采样点
		double uMin, uMax;
		Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, uMin, uMax);
		if (curve.IsNull()) {
			return false; // 曲线为空
		}
		double step = (uMax - uMin) / sampleCount;
		double u = uMin;
		// 不采样edge的终点
		for (int j = 0; j < sampleCount; j++) {
			gp_Pnt pnt;
			curve->D0(u, pnt);
			samplePnts.push_back(pnt);
			u += step;
		}
	}

	if (samplePnts.empty()) {
		return false; // 没有采样点
	}

	samplePnts.push_back(BRep_Tool::Pnt(TopExp::LastVertex(continueEdge->edges.back(), true)));
	return true;
}

TopoDS_Shape CuttingRailsImpl::ShowEdge(const Handle(Geom_Curve)& curve, double first, double last) {
	BRepBuilderAPI_MakeEdge makeEdge(curve, first, last);
	return makeEdge.Edge();
}

bool CuttingRailsImpl::HasVRotationalAxis(const Handle(Geom_Surface)& surface, double uMin, double uMax, double vMin, double vMax, gp_Vec& direction) {
	// 检查多个v参数下的u等参线
	int numSections = 6;
	std::vector<gp_Vec> normals;
	for (int i = 0; i < numSections; i++) {
		double v = vMin + (vMax - vMin) * i / (numSections - 1);
		Handle(Geom_Curve) curve = surface->VIso(v);
		//double f = curve->FirstParameter();
		//double l = curve->LastParameter();
		//TopoDS_Shape isoEdge = ShowEdge(curve, std::max(uMin, f), std::min(uMax, l));
		gp_Pnt center;
		gp_Vec normal;
		if (IsPlanerCurve(curve, uMin, uMax, normal) > 0) {
			normals.push_back(normal);
		}		
	}
	if (normals.size() < 3) {
		return false;
	}
	// 检查圆心是否共线
	direction = normals[0];
	for (size_t i = 1; i < normals.size(); i++) {
		if (!direction.IsParallel(normals[i], M_PI / 36))
			return false;
	}

	return true;
}

bool CuttingRailsImpl::HasURotationalAxis(const Handle(Geom_Surface)& surface, double uMin, double uMax, double vMin, double vMax, gp_Vec& direction) {
	// 检查多个u参数下的v等参线
	int numSections = 6;
	std::vector<gp_Vec> normals;
	for (int i = 0; i < numSections; i++) {
		double u = uMin + (uMax - uMin) * i / (numSections - 1);
		Handle(Geom_Curve) curve = surface->UIso(u);
		//double f = curve->FirstParameter();
		//double l = curve->LastParameter();
		//TopoDS_Shape isoEdge = ShowEdge(curve, std::max(vMin, f), std::min(vMax, l));
		gp_Vec normal;
		if (IsPlanerCurve(curve, vMin, vMax, normal) > 0) {
			normals.push_back(normal);
		}
	}

	if (normals.size() < 3) {
		return false;
	}

	direction = normals[0];
	for (size_t i = 1; i < normals.size(); i++) {
		if (!direction.IsParallel(normals[i], M_PI / 36))
			return false;
	}

	return true;
}

bool CuttingRailsImpl::HasSphereNearExtrudeDir(const TopoDS_Face& face, const gp_Vec& sphereNormal, gp_Vec& dir, double tol, double squareTol) {
	// 找出最长的连续边, 并对其进行判断是否为圆弧
	TopoDS_Shape mergedFace = MergeSameVerticesForFaces(face);
	TopTools_IndexedDataMapOfShapeShape faceMapping;
	mergedFace = FixFaceWiresOrder(mergedFace, faceMapping);
	std::vector<ContinueEdgePtr> allWireContinuesEdges;
	
	// std::cout << "all continues count: " << allWireContinuesEdges.size() << std::endl;
	std::vector<ContinueEdgePtr> firstUnmatched;
	std::vector<ContinueEdgePtr> unmatchedAfterSeamEdge;
	std::vector<ContinueEdgePtr> unmatchedAfterMergeEdge;

	TopTools_IndexedMapOfShape softVertices;
	EdgeSplitMap edgeSplitParams;
	GetSoftVertices(mergedFace, tol, softVertices, edgeSplitParams);
	mergedFace = SplitShapeEdges(mergedFace, edgeSplitParams);
	//TopoDS_Shape sofeVertices = showTopos(softVertices);
	std::shared_ptr<KDTree> tree = MakeKDTree(mergedFace);
	MakeContinueEdges(tree, mergedFace, allWireContinuesEdges, softVertices, tol);
	double maxLen = allWireContinuesEdges.front()->len;
	int maxLenIndex = 0;
	for (int i = 1, iLen = allWireContinuesEdges.size(); i < iLen; i++) {
		double curLen = allWireContinuesEdges[i]->len;
		if (curLen > maxLen) {
			maxLen = curLen;
			maxLenIndex = i;
		}
	}

	ContinueEdgePtr maxContinueEdge = allWireContinuesEdges[maxLenIndex];
	std::vector<gp_Pnt> samplePnts;
	if (!SampleContinue(maxContinueEdge, 100, samplePnts)) {
		return false;
	}

	if (samplePnts.size() < 3) {
		return false; // 采样点不足
	}

	//gp_Pnt center;
	//if (!IsCircularCurve(samplePnts, center, dir, squareTol)) {
	//	return false; // 不是圆弧
	//}
	if (IsPlanerCurve(samplePnts, dir, squareTol) < 1) {
		// 计算每个点投影到normal的参数
		std::vector<std::pair<double, int>> dotValues;
		gp_Vec vec;
		for (int i = 0, iLen = samplePnts.size(); i < iLen; i++) {
			vec.SetXYZ(samplePnts[i].XYZ());
			double dotValue = vec.Dot(sphereNormal);
			dotValues.push_back({ dotValue, i });
		}

		std::sort(dotValues.begin(), dotValues.end(),
			[](const std::pair<double, int>& a, const std::pair<double, int>& b) {
				return a.first < b.first;
			});

		if (IsPlanerCurve(
			{ samplePnts[dotValues.front().second],
			samplePnts[dotValues[dotValues.size() / 2].second],
			samplePnts[dotValues.back().second] }, dir, squareTol) < 1) {
			return false; // 不是平面曲线
		}
	}
	return true;
}

bool CuttingRailsImpl::GetSpecificSurfExtrudeDir(const Handle(Geom_Surface)& surf, gp_Vec& dir)
{
	// 判断是否为旋转面
	if (surf->IsKind(STANDARD_TYPE(Geom_SurfaceOfRevolution))) {
		// 获取旋转轴
		Handle(Geom_SurfaceOfRevolution) revSurf =
			Handle(Geom_SurfaceOfRevolution)::DownCast(surf);
		gp_Ax1 axis = revSurf->Axis();
		dir = axis.Direction();
		return true;
	}

	// 判断是否为圆柱面
	if (surf->IsKind(STANDARD_TYPE(Geom_CylindricalSurface))) {
		Handle(Geom_CylindricalSurface) cylSurf =
			Handle(Geom_CylindricalSurface)::DownCast(surf);
		gp_Ax1 axis = cylSurf->Axis();
		dir = axis.Direction();
		return true;
	}

	// 判断是否为圆锥面
	if (surf->IsKind(STANDARD_TYPE(Geom_ConicalSurface))) {
		Handle(Geom_ConicalSurface) conSurf =
			Handle(Geom_ConicalSurface)::DownCast(surf);
		gp_Ax1 axis = conSurf->Axis();
		dir = axis.Direction();
		return true;
	}

	// 判断是否为拉伸面
	if (surf->IsKind(STANDARD_TYPE(Geom_SurfaceOfLinearExtrusion))) {
		Handle(Geom_SurfaceOfLinearExtrusion) extSurf =
			Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(surf);
		// 拉伸面只有在基准曲线为圆时才有旋转轴
		Handle(Geom_Curve) baseCurve = extSurf->BasisCurve();
		if (baseCurve->IsKind(STANDARD_TYPE(Geom_Circle))) {
			gp_Dir extDir = extSurf->Direction();
			Handle(Geom_Circle) circle = Handle(Geom_Circle)::DownCast(baseCurve);
			gp_Ax1 circleAxis = circle->Axis();
			// 判断拉伸方向是否与圆的轴线平行
			if (extDir.IsParallel(circleAxis.Direction(), Precision::Angular())) {
				return true;
			}
		}
	}
	return false;
}

// 辅助函数：解 (A - lambda * I) * x = 0，并把结果存到 vecOut
// 若出现严重数值问题(如pivot很小)，返回false
bool CuttingRailsImpl::solveEigenVector3x3(const std::array<double, 9>& A,
	double lambda,
	std::array<double, 3>& vecOut)
{
	// 构造 B = A - λI
	std::array<double, 9> B = A;
	B[0] -= lambda;  // A[0,0]
	B[4] -= lambda;  // A[1,1]
	B[8] -= lambda;  // A[2,2]

	// 3x3简易高斯消去(带行pivot)，若出现数值奇异则返回false
	const double TINY = 1e-14;

	// 第1行pivot
	int pivotRow = 0;
	double maxVal = std::fabs(B[0]);
	for (int r = 1; r < 3; r++) {
		double tmp = std::fabs(B[r * 3 + 0]);
		if (tmp > maxVal) {
			pivotRow = r;
			maxVal = tmp;
		}
	}
	// 如果pivot极小，说明矩阵高度奇异 => 返回false
	if (maxVal < TINY) return false;

	// 若pivotRow!=0则交换
	if (pivotRow != 0) {
		for (int c = 0; c < 3; c++) {
			std::swap(B[c], B[pivotRow * 3 + c]);
		}
	}
	double pivotVal = B[0];
	// 用第1行消去第2,3行的第一列
	for (int r = 1; r < 3; r++) {
		double factor = B[r * 3 + 0] / pivotVal;
		for (int c = 0; c < 3; c++) {
			B[r * 3 + c] -= factor * B[c];
		}
	}

	// 第2行pivot
	pivotRow = 1;
	maxVal = std::fabs(B[1 * 3 + 1]);
	if (std::fabs(B[2 * 3 + 1]) > maxVal) {
		pivotRow = 2;
		maxVal = std::fabs(B[2 * 3 + 1]);
	}
	if (maxVal < TINY) return false;

	// 若pivotRow!=1则交换第1、2行
	if (pivotRow != 1) {
		for (int c = 0; c < 3; c++) {
			std::swap(B[1 * 3 + c], B[2 * 3 + c]);
		}
	}

	double pivotVal2 = B[1 * 3 + 1];
	double factor2 = B[2 * 3 + 1] / pivotVal2;
	for (int c = 1; c < 3; c++) {
		B[2 * 3 + c] -= factor2 * B[1 * 3 + c];
	}

	// 回代：假定 x2=1 然后算 x1,x0
	// 如果 B[2,2] 也接近0，则说明有更高的自由度，这里可直接设 x2=1
	// 只要不全为0即可
	double pivot22 = std::fabs(B[2 * 3 + 2]);
	double x2 = (pivot22 < TINY) ? 1.0 : 1.0; // 反正设成1就行
	double x1 = 0.0;
	{
		double diag11 = B[1 * 3 + 1];
		double b12 = B[1 * 3 + 2];
		if (std::fabs(diag11) > TINY) {
			x1 = -b12 * x2 / diag11;
		}
		else {
			// 如果这里也很小，说明 x1自由，可随意设0
			x1 = 0.0;
		}
	}
	double x0 = 0.0;
	{
		double diag00 = B[0];
		double b01 = B[1];
		double b02 = B[2];
		if (std::fabs(diag00) > TINY) {
			x0 = -(b01 * x1 + b02 * x2) / diag00;
		}
		else {
			x0 = 0.0;
		}
	}

	// 归一化
	double normVal = std::sqrt(x0 * x0 + x1 * x1 + x2 * x2);
	if (normVal < TINY) return false;

	vecOut = { x0 / normVal, x1 / normVal, x2 / normVal };
	return true;
}

// 3x3特征值分解：
// - 若成功，返回true，并在 outEigenPairs 中填入 3 个条目(对复根仅填实部虚部)
// - 若出现数值问题(如奇异，NaN等)，返回false
bool CuttingRailsImpl::EigenDecompose3x3(const std::array<double, 9>& A,
	std::vector<EigenPair>& outEigenPairs)
{
	// 提取元素
	double a00 = A[0], a01 = A[1], a02 = A[2],
		a10 = A[3], a11 = A[4], a12 = A[5],
		a20 = A[6], a21 = A[7], a22 = A[8];

	// 特征多项式 λ^3 + c2 λ^2 + c1 λ + c0 = 0
	double c2 = -(a00 + a11 + a22);
	double c1 = (a00 * a11 + a00 * a22 + a11 * a22)
		- (a01 * a10 + a02 * a20 + a12 * a21);
	double detA =
		a00 * (a11 * a22 - a12 * a21)
		- a01 * (a10 * a22 - a12 * a20)
		+ a02 * (a10 * a21 - a11 * a20);
	double c0 = -detA;

	double inv3 = 1.0 / 3.0;
	double c2_3 = c2 * inv3;
	double p = c1 - 3.0 * c2_3 * c2_3;
	double q = 2.0 * std::pow(c2_3, 3) - c2_3 * c1 + c0;
	double disc = 0.25 * q * q + (1.0 / 27.0) * p * p * p;

	// 若数值过大或NaN
	if (std::isnan(disc) || std::fabs(disc) > 1e20) {
		return false;
	}

	outEigenPairs.clear();

	const double EPS = 1e-14;
	if (disc > EPS) {
		// 1个实根 + 1对复共轭
		double sqd = std::sqrt(disc);
		double alpha = -0.5 * q + sqd;
		double beta = -0.5 * q - sqd;

		auto cbrt_sign = [&](double val) {
			return val >= 0.0 ? std::cbrt(val) : -std::cbrt(-val);
			};
		double A_ = cbrt_sign(alpha);
		double B_ = cbrt_sign(beta);

		double realRoot = A_ + B_;
		double real2 = -0.5 * (A_ + B_);
		double imag2 = (std::sqrt(3.0) / 2.0) * (A_ - B_);

		realRoot -= c2_3;
		real2 -= c2_3;

		// 实根
		EigenPair epR;
		epR.isComplex = false;
		epR.realPart = realRoot;
		outEigenPairs.push_back(epR);

		// 复根
		EigenPair epC1, epC2;
		epC1.isComplex = true;
		epC1.realPart = real2;
		epC1.imagPart = imag2;

		epC2.isComplex = true;
		epC2.realPart = real2;
		epC2.imagPart = -imag2;

		outEigenPairs.push_back(epC1);
		outEigenPairs.push_back(epC2);
	}
	else {
		// 3个实根
		double rho = std::sqrt(-(1.0 / 27.0) * p * p * p);
		if (rho < EPS) {
			// 极端情况 => 返回false
			return false;
		}
		double theta = std::acos(-0.5 * q / rho);
		double two_r = 2.0 * std::cbrt(rho);

		double x1 = two_r * std::cos(theta / 3.0);
		double x2 = two_r * std::cos((theta + 2.0 * M_PI) / 3.0);
		double x3 = two_r * std::cos((theta + 4.0 * M_PI) / 3.0);

		x1 -= c2_3;
		x2 -= c2_3;
		x3 -= c2_3;

		EigenPair ep1, ep2, ep3;
		ep1.isComplex = false; ep1.realPart = x1;
		ep2.isComplex = false; ep2.realPart = x2;
		ep3.isComplex = false; ep3.realPart = x3;
		outEigenPairs.push_back(ep1);
		outEigenPairs.push_back(ep2);
		outEigenPairs.push_back(ep3);
	}

	// 若有实特征值，则解向量
	for (auto& ep : outEigenPairs) {
		if (!ep.isComplex) {
			std::array<double, 3> vec = { 0,0,0 };
			bool ok = solveEigenVector3x3(A, ep.realPart, vec);
			if (!ok) {
				// 若solve失败，直接返回false
				return false;
			}
			ep.eigenVector = vec;
		}
	}
	return true;
}

bool CuttingRailsImpl::GetFaceExtrudeDirByPCA(const TopoDS_Face& face, gp_Vec& dir) {
	// 使用 BRepBuilderAPI_Copy 对 face 进行复制
	BRepBuilderAPI_Copy copyTool(face);
	TopoDS_Face copiedFace = TopoDS::Face(copyTool.Shape());
	TopLoc_Location theLocation;
	Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(copiedFace, theLocation);
	if (triangulation.IsNull()) {
		BRepMesh_IncrementalMesh mesh(copiedFace, 0.05);
		triangulation = BRep_Tool::Triangulation(copiedFace, theLocation);
	}

	if (triangulation.IsNull()) {
		return false; // 无法获取三角化数据
	}

	// 获取三角化数据的顶点
	std::vector<gp_Pnt> points;
	double centerX = 0;
	double centerY = 0;
	double centerZ = 0;

	for (int i = 1; i <= triangulation->NbNodes(); ++i) {
		gp_Pnt p = triangulation->Node(i);
		points.push_back(p);
		centerX += p.X();
		centerY += p.Y();
		centerZ += p.Z();
	}

	int numPoints = triangulation->NbNodes();
	if (numPoints == 0) {
		return false; // 无顶点
	}

	centerX /= numPoints;
	centerY /= numPoints;
	centerZ /= numPoints;

	std::array<double, 9> covMat = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	for (auto& p : points) {
		// 当前点相对于质心的偏移向量
		double dx = p.X() - centerX;
		double dy = p.Y() - centerY;
		double dz = p.Z() - centerZ;
		// 累加到 covMat
		covMat[0] += dx * dx;
		covMat[1] += dx * dy;
		covMat[2] += dx * dz;

		covMat[3] += dy * dx;
		covMat[4] += dy * dy;
		covMat[5] += dy * dz;

		covMat[6] += dz * dx;
		covMat[7] += dz * dy;
		covMat[8] += dz * dz;
	}

	for (int i = 0; i < 9; ++i) {
		covMat[i] /= numPoints;
	}

	// 计算协方差矩阵的特征值和特征向量
	std::vector<EigenPair> eigenPairs;
	if (EigenDecompose3x3(covMat, eigenPairs)) {
		// 找到最大的特征值对应的特征向量
		EigenPair maxEigen = eigenPairs[0];
		for (const auto& ep : eigenPairs) {
			if (ep.realPart > maxEigen.realPart) {
				maxEigen = ep;
			}
		}
		if (maxEigen.IsValid()) {
			dir = gp_Vec(maxEigen.eigenVector[0], maxEigen.eigenVector[1], maxEigen.eigenVector[2]);
			dir.Normalize();
			// 设置 extrude 方向
			return true; // 成功获取 extrude 方向
		}
	}

	return false;
}

bool CuttingRailsImpl::GetExtrudeDirByPtsPCA(const TopoDS_Shape& shape, gp_Vec& dir) {
	// 使用 BRepBuilderAPI_Copy 对 face 进行复制
	TopTools_IndexedMapOfShape vertices;
	TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);

	// 获取三角化数据的顶点
	std::vector<gp_Pnt> points;
	points.reserve(vertices.Size());
	double centerX = 0;
	double centerY = 0;
	double centerZ = 0;

	for (int i = 1; i <= vertices.Size(); i++) {
		points.push_back(BRep_Tool::Pnt(TopoDS::Vertex(vertices.FindKey(i))));
	}

	int numPoints = points.size();
	if (numPoints == 0) {
		return false; // 无顶点
	}

	centerX /= numPoints;
	centerY /= numPoints;
	centerZ /= numPoints;

	std::array<double, 9> covMat = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	for (auto& p : points) {
		// 当前点相对于质心的偏移向量
		double dx = p.X() - centerX;
		double dy = p.Y() - centerY;
		double dz = p.Z() - centerZ;
		// 累加到 covMat
		covMat[0] += dx * dx;
		covMat[1] += dx * dy;
		covMat[2] += dx * dz;

		covMat[3] += dy * dx;
		covMat[4] += dy * dy;
		covMat[5] += dy * dz;

		covMat[6] += dz * dx;
		covMat[7] += dz * dy;
		covMat[8] += dz * dz;
	}

	for (int i = 0; i < 9; ++i) {
		covMat[i] /= numPoints;
	}

	// 计算协方差矩阵的特征值和特征向量
	std::vector<EigenPair> eigenPairs;
	if (EigenDecompose3x3(covMat, eigenPairs)) {
		// 找到最大的特征值对应的特征向量
		EigenPair maxEigen = eigenPairs[0];
		for (const auto& ep : eigenPairs) {
			if (ep.realPart > maxEigen.realPart) {
				maxEigen = ep;
			}
		}
		if (maxEigen.IsValid()) {
			dir = gp_Vec(maxEigen.eigenVector[0], maxEigen.eigenVector[1], maxEigen.eigenVector[2]);
			dir.Normalize();
			// 设置 extrude 方向
			return true; // 成功获取 extrude 方向
		}
	}

	return false;
}

bool CuttingRailsImpl::GetNearExtrudeDirection(const TopoDS_Shape& sortedShape, const std::vector<double>& sortedFaceAreas, gp_Vec& dir, double tol) {
	TopExp_Explorer exp(sortedShape, TopAbs_FACE);
	bool find = false;
	double squareTol = 1.0e-4;
	gp_Vec sphereNormal;
	// 选择面积最大的5个面, 优先判断柱面、拉伸面、锥面、旋转面
	std::vector<TopoDS_Face> faces;
	std::vector<Handle(Geom_Surface)> surfs;
	bool checkSpecificFirst = true;
	for (int i = 0; i < 5 && exp.More(); exp.Next(), i++) {
		TopoDS_Face face = TopoDS::Face(exp.Current());
		Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
		if (surf->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface))) {
			surf = Handle(Geom_RectangularTrimmedSurface)::DownCast(surf)->BasisSurface();
		}

		faces.push_back(face);
		surfs.push_back(surf);

		if (i > 0 && sortedFaceAreas[i] < sortedFaceAreas[i - 1] * 0.5) {
			checkSpecificFirst = false;
			continue;
		}
		
		if (checkSpecificFirst && GetSpecificSurfExtrudeDir(surf, dir)) {
			return true;
		}
	}

	for (int i = 0; i < faces.size(); i++) {
		TopoDS_Face face = faces[i];
		
		Handle(Geom_Surface) surf = surfs[i];

		if (GetSpecificSurfExtrudeDir(surf, dir)) {
			return true;
		}

		// 判断是否为球面
		bool isSphere = false;
		if (surf->IsKind(STANDARD_TYPE(Geom_SphericalSurface))) {
			// 球面有无数个旋转轴,都通过球心
			isSphere = true;
			sphereNormal = Handle(Geom_SphericalSurface)::DownCast(surf)->Axis().Direction();
		}

		if (!isSphere) {
			gp_Vec dir1, dir2;
			double umin, umax, vmin, vmax;
			BRepTools::UVBounds(face, umin, umax, vmin, vmax);
			bool hasUDir = HasURotationalAxis(surf, umin, umax, vmin, vmax, dir1);
			bool hasVDir = HasVRotationalAxis(surf, umin, umax, vmin, vmax, dir2);
			if (hasUDir && hasVDir) {
				gp_Vec norm = dir1.Crossed(dir2);
				if (norm.SquareMagnitude() > gp::Resolution()) {
					sphereNormal = norm.Normalized();
					isSphere = true;
				}
			}
			else if (hasUDir) {
				dir = dir1;
				return true;
			}
			else if (hasVDir) {
				dir = dir2;
				return true;
			}
		}

		if (isSphere && HasSphereNearExtrudeDir(face, sphereNormal, dir, tol, squareTol)) {
			// 法向修正
			// 计算与三个轴的夹角（实际比较的是 cos 值的绝对值，越大越接近）
			Standard_Real dotX = std::abs(dir.Dot(gp::DX()));
			Standard_Real dotY = std::abs(dir.Dot(gp::DY()));
			Standard_Real dotZ = std::abs(dir.Dot(gp::DZ()));

			// 找出最大分量所对应的轴方向
			if (dotX >= dotY && dotX >= dotZ) {
				dir = gp_Vec(gp::DX());
			}
			else if (dotY >= dotX && dotY >= dotZ) {
				dir = gp_Vec(gp::DY());
			}
			else {
				dir = gp_Vec(gp::DZ());
			}
			return true;
		}
	}

	//// 当没有找到一个合适的旋转轴时，尝试使用第一个面的离散网格做PCA分析, 作为拉伸方向, 失败返回false
	//find = GetFaceeExtrudeDirByPCA(maxFace, dir);
	return find;
}
