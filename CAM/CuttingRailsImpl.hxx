
#ifndef _CUTTINGRAILSIMPL_HeaderFile
#define _CUTTINGRAILSIMPL_HeaderFile

#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <Geom_Surface.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Compound.hxx>
#include <Geom_Curve.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_IndexedDataMapOfShapeShape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <BRepClass3d_SolidExplorer.hxx>
#include <TColStd_ListOfReal.hxx>

#include <array>
#include <vector>
#include <set>

#include "../GEOM/Geom_CurveOnSurface.hxx"
#include "../UTILS/KDtree.hxx"
#include "./CuttingRailsOption.hxx"

class CuttingRailsImpl {
public:
	bool hasOnVertexOnOnEdge = false;
	enum AdjecentType {
		OUTER,
		INNER,
		BOTH,
		UNKNOWNAdjecentType
	};

	enum EdgeType {
		TU,
		AO,
		SMOOTH,
		NOTKNOWN
	};

	enum ErrorCode {
		ThreeCoedge,
	};

	// 用 shape.IsSame 判定是否“相同”
	struct ShapeEqual {
		bool operator()(const TopoDS_Shape& a, const TopoDS_Shape& b) const {
			return a.IsSame(b);
		}
	};

	// 用于存储每一段连续区间的信息
	struct RunInfo {
		int value;  // 数字
		int start;  // 该段连续区间的起始下标
		int end;    // 该段连续区间的结束下标
	};

	struct PointWithParam {
	public:
		gp_Pnt point;
		double t;

		PointWithParam(double param, const gp_Pnt& p) : point(p), t(param) {}
	};

	// 用于返回特征值和向量
	struct EigenPair {
		bool isComplex = false;
		double realPart = 0.0;
		double imagPart = 0.0;
		// 对于实特征值，这里存特征向量；对复特征值则不使用
		std::array<double, 3> eigenVector = { 0,0,0 };

		bool IsValid() {
			return eigenVector[0] * eigenVector[0] + eigenVector[1] * eigenVector[1] + eigenVector[2] * eigenVector[2] > 1.0e-12;
		}
	};

	struct ContinueEdge {
		std::vector<TopoDS_Edge> edges;
		gp_Pnt start;
		gp_Pnt end;
		gp_Pnt mid;
		gp_Vec midVec;
		gp_Dir midNormal;
		TopoDS_Shape face;
		TopoDS_Shape twinface;
		bool bOuterLoop;
		double len;
		EdgeType edgetType;
		bool degenerated;
		int faceIndex = -1;
		std::vector<int> edgesIndex;
		bool bClosed;
		TopoDS_Compound compedges;
		std::vector<double> increaseLen;
		// 用于merge 构造
		ContinueEdge(const std::vector<TopoDS_Edge>& es, double length);

		ContinueEdge(const std::vector<TopoDS_Edge>& es, const TopoDS_Shape& ff, bool bOutLoop, double length, int index);

		ContinueEdge(const std::vector<TopoDS_Edge>& es, const TopoDS_Shape& f, bool bOutLoop, int faceIndex);

	public:

		std::shared_ptr<ContinueEdge> Reversed();

		TopoDS_Shape GetCompound();

		double GetLength(std::vector<double>& increase);

		void GetNormalOfMiddleEdge(const int midIndex);

		gp_Pnt2d uvAtStart();

		bool NormalAt(const gp_Pnt2d& pt2d, gp_Vec& normal);

		gp_Pnt2d uvAtEnd();

		gp_Vec dirAtStart();

		gp_Vec dirAtEnd();

		TopoDS_Vertex GetStartVertex();

		TopoDS_Vertex GetEndVertex();

		bool UpdateMidPt(const gp_Vec& midTangent, const gp_Pnt2d& midUV);
	};

	typedef std::shared_ptr<ContinueEdge> ContinueEdgePtr;

	typedef std::vector<std::pair<gp_Pnt, int>> VertexToContinueIndexMap;

	typedef std::vector<std::pair<std::vector<ContinueEdgePtr>, bool>> CoedgeInfo;

	typedef NCollection_IndexedDataMap
		<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> shape2shapeMap;

	typedef std::vector <std::pair<ContinueEdgePtr, ContinueEdgePtr>> EdgeVector;

	typedef std::shared_ptr<EdgeVector> EdgeVectorPtr;

	typedef NCollection_IndexedDataMap
		<TopoDS_Shape, EdgeVectorPtr, TopTools_ShapeMapHasher> Face2CoedgesMap;

	typedef NCollection_IndexedDataMap
		<TopoDS_Shape, std::vector<size_t>, TopTools_ShapeMapHasher> VertexToContinueMap;

	typedef NCollection_IndexedDataMap
		<TopoDS_Shape, std::vector<ContinueEdgePtr>, TopTools_ShapeMapHasher> FaceToContinuesMap;

	typedef NCollection_IndexedDataMap<TopoDS_Edge, std::vector<double>> EdgeSplitMap;

	typedef NCollection_IndexedDataMap<TopoDS_Shape, int> FaceIndexMap;

public:
	// 测试可视化
	const char* BRepTools_Write(const char* theFileStr, void* theShapePtr);

	TopoDS_Shape showTopos(const TopTools_IndexedMapOfShape& shapes);

	void PrintFaceEdge(const std::vector<ContinueEdgePtr>& wireContinueEdges, int fIndex = 0);

	TopoDS_Shape ShowEdges(const TopTools_MapOfShape& edges);

	TopoDS_Shape ShowEdges(const TopTools_IndexedMapOfShape& edges);

	TopoDS_Compound showEdgesType(const CoedgeInfo& coedgeInfoVector);

	TopoDS_Shape showAllContinuesEdges(const std::vector<ContinueEdgePtr>& allWireContinuesEdges);

	TopoDS_Shape ShowEdge(const Handle(Geom_Curve)& curve, double first, double last);

public:
	bool HasFace(const TopoDS_Shape& shape);

	TopoDS_Compound copyShape(const TopoDS_Shape& shape);

	void GetCuttingRails(const TopoDS_Shape& shape, TopoDS_Shape& outerBoundFaces, TopoDS_Shape& innerBoundFaces, TopoDS_Compound& unMatchedEdges,
		std::vector<std::vector<ContinueEdgePtr>>& loops, std::vector<std::vector<ContinueEdgePtr>>& twinLoops, const CuttingRailsOption& option, double tol);

	TopoDS_Shape SortFaceByArea(const TopoDS_Shape& shape, std::vector<double>& areas);

	//要求面积排过序
	bool GetExtrudeDirection(const TopoDS_Shape& shape, gp_Vec& dir);

	bool GetNearExtrudeDirection(const TopoDS_Shape& sortedShape, const std::vector<double>& sortedFaceAreas, gp_Vec& dir, double tol = 0.1);

	bool FindDirectionBySurfaceNormalCross(const TopoDS_Shape& shape, gp_Vec& direction);

	bool GetSortedShapeExtrudeDir(const TopoDS_Shape& shape, const std::vector<double>& sortedFaceAreas, gp_Vec& dir);

private:

	TopoDS_Shape GetFirstSubShape(const TopoDS_Shape& shape);

	Standard_Boolean FaceNormal(const TopoDS_Face& aF,
		const Standard_Real U,
		const Standard_Real V,
		gp_Dir& aDN);

	double dist2VecVec(const point_t& p1, const point_t& pt);

	std::shared_ptr<KDTree> MakeKDTree(const TopoDS_Shape& shape);

	bool IsPtBNearestPtA(const std::shared_ptr<KDTree>& tree, const point_t& ptA, const point_t& ptB, double tol);

	gp_Vec DirAtStart(const TopoDS_Edge& edge);

	bool IsSmoothAtVertex(const TopTools_IndexedDataMapOfShapeListOfShape& VEMap, const TopoDS_Vertex& vertex);

	bool MakeContinueEdges(const std::shared_ptr<KDTree>& tree, const TopoDS_Shape& shape,
		std::vector<ContinueEdgePtr>& wireContinueEdges, TopTools_IndexedMapOfShape& softVs, double tol);

	void GetEdgeCurveTangent(const TopoDS_Shape& edge, const Handle(Geom_Curve)& curve, double param, gp_Vec& tangent);

	double PointDistanceToCompCurve(const gp_Pnt& p, const ContinueEdgePtr& compCurve, gp_Vec& refMidVec, gp_Pnt2d& refMidUV);

	int NbChildren(const TopoDS_Shape& shape);

	TopoDS_Compound GetUnMatchedEdges(const std::vector<ContinueEdgePtr>& CCUnMatched, const double tol);

	int MatchBFromA(const std::shared_ptr<KDTree>& tree,
		const VertexToContinueMap& vertexToContinueMap,
		const size_t ce1Index,
		const std::vector<ContinueEdgePtr>& continueEdges,
		double tol,
		double midTol, bool& sameDir, gp_Vec& matchedMidVec, gp_Pnt2d& matchedMidUV);

	bool MatchCoedges(const std::shared_ptr<KDTree>& tree,
		const std::vector<ContinueEdgePtr>& continueEdges,
		CoedgeInfo& CCMatchedMap,
		std::vector<ContinueEdgePtr>& CCUnMatched,
		double tol,
		TopoDS_Compound& unMatchedEdges);

	std::vector<ContinueEdgePtr> GetCandidateRails(const TopoDS_Shape& outBoundFaces, const Face2CoedgesMap& faceCoedgesMap, bool needCuttingFace);

	bool Contains(const VertexToContinueIndexMap& vertexMap, const gp_Pnt& point, const double squareTol);

	std::vector<int> GetValue(const VertexToContinueIndexMap& vertexMap, const gp_Pnt& point, const double squareTol);

	VertexToContinueIndexMap BuildVertexToContinueIndexMap(const std::vector<ContinueEdgePtr>& continues);
	// 查找下一条边
	int FindNextEdge(const gp_Pnt& point,
		const VertexToContinueIndexMap& vertexMap,
		const std::set<int>& usedEdges,
		double squareTol);

	std::vector<int> SearchLoop(const int startIndex,
		const VertexToContinueIndexMap& vertexMap,
		const std::vector<ContinueEdgePtr>& toBeSearched,
		std::set<int>& usedEdges,
		const TopTools_IndexedMapOfShape& borderVertices,
		double squareTol);

	void GetBorderVertex(const std::vector<ContinueEdgePtr>& continues, TopTools_IndexedMapOfShape& borderVertices);

	std::vector<ContinueEdgePtr> FindLoops(const std::vector<ContinueEdgePtr>& continues, std::vector<std::vector<ContinueEdgePtr>>& loops, const double tolerance);

	TopoDS_Shape FixFaceWiresOrder(const TopoDS_Shape& shape, TopTools_IndexedDataMapOfShapeShape& faceMapping);

	// 平面返回normal, 拉伸面返回 dir
	bool IsExtrudeFace(const TopoDS_Face& face, gp_Vec& dir1, gp_Vec& dir2, double theAngularTolerance);

	TopoDS_Shape SortShellByArea(const TopoDS_Shape& shape);

	TopoDS_Shape FixFaceOritation(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfos, TopTools_IndexedMapOfShape& reversedFaces);

	void ClassifyInOutGroups(const TopoDS_Compound& groups,
		TopoDS_Shape& outerBoundFaces, TopoDS_Shape& innerBoundFaces, const bool hasExtrudeDir, const gp_Vec& extrudeDir, const SingleBorderType& singleBorderType);

	bool ClassifyInOutGroups(const TopoDS_Shape& groups, TopoDS_Shape& outerBoundFaces, TopoDS_Shape& innerBoundFaces);

	TopoDS_Shape GetBiggestGroup(const TopoDS_Shape& groups);

	Standard_Real ComputeThicknessByRay(const TopoDS_Shape& shape, const TopoDS_Face& skipFace, const CoedgeInfo& coedgeInfoVector);

	bool IsLeftRightOnSameSide(const EdgeVectorPtr& coedges, const ContinueEdgePtr& coedge, double tol);

	bool IsSameClassification(const std::pair<ContinueEdgePtr, ContinueEdgePtr>& edge, const Face2CoedgesMap& faceCoedgesMap, double thinDist, double tol);

	bool IsSameClassificationWithoutLength(const std::pair<ContinueEdgePtr, ContinueEdgePtr>& edge, const Face2CoedgesMap& faceCoedgesMap, double tol);

	void AddFacesToGroup(const TopoDS_Shape& shape, TopTools_IndexedMapOfShape& faceGroup);

	TopoDS_Face GetOtherSideFace(const TopoDS_Shape& sommothFace, BRepClass3d_SolidExplorer& aSE, double& minDist);

	bool IsInnerOuterShell(const TopoDS_Shell& shell, const Face2CoedgesMap& faceCoedgesMap, double thinDist, double tol);

	bool HasFaceInShell(const TopoDS_Face& targetFace, const TopoDS_Shell& shell);

	TopoDS_Shell GetRelationPatchs(const TopoDS_Shape& otherSideFace, const TopoDS_Shape& fixedShape);

	void GetAdjecentFaces(const TopoDS_Shape& shell, const Face2CoedgesMap& face2Coedges,
		TopTools_IndexedMapOfShape& adjecentFaces);

	void GetAdjecentSmoothFaces(const TopoDS_Shape& shell, const Face2CoedgesMap& face2Coedges,
		TopTools_IndexedMapOfShape& adjecentSmoothFaces, double tol);

	bool HasCommonFace(const TopoDS_Shape& patchShell, const TopTools_IndexedMapOfShape& adjecentFaces);

	void GetSameClassificationFromFuzzyFaces(const TopoDS_Face& face, const TopTools_IndexedMapOfShape& fuzzyFaces, const Face2CoedgesMap& faceCoedgesMap, TopTools_IndexedMapOfShape& faceGroup, double thinDist, double tol);

	TopoDS_Compound MakeSameClassificationForFuzzyFaces(const TopTools_IndexedMapOfShape& fuzzyFaces, const Face2CoedgesMap& faceCoedgesMap, double thinDist, double tol);

	TopoDS_Shell GetPatchFromFace(const TopoDS_Compound& comp, const TopoDS_Face& face);

	void GetInnerOuterShellByRatio(const TopoDS_Shape& fixedShape, const Face2CoedgesMap& faceCoedgesMap,
		TopTools_IndexedMapOfShape& innerouterFaces, TopTools_IndexedMapOfShape& cuttingFaces,
		TopTools_IndexedMapOfShape& fuzzyFaces,
		std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair,
		bool bExtrudeShape, double thinDist, double tol);

	void GetFaceContinuesEdges(const TopoDS_Shape& face, const CoedgeInfo& coedgeInfos, const EdgeVectorPtr& faceCoedges);

	void MakeEdgeType(const CoedgeInfo& coedges);

	// 在shape里找相邻
	void GetSameClassificationPatches(const TopoDS_Shape& face, const TopTools_IndexedMapOfShape& allFaces, const Face2CoedgesMap& faceCoedgesMap,
		TopTools_IndexedMapOfShape& faceGroup, double thinDist, double tol);

	TopoDS_Shape MakeSameClassificationPatchs(const TopTools_IndexedMapOfShape& shape, const Face2CoedgesMap& faceCoedgesMap,
		double thinDist, double tol, bool bConsiderLength);

	bool IsExclusivePair(const TopoDS_Face& face1, const TopoDS_Face& face2,
		const std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair);

	TopoDS_Shell MakeFaceGroup(const TopoDS_Face& face, const Face2CoedgesMap& faceCoedgesMap,
		TopTools_IndexedMapOfShape& visitedFaceGroup,
		const TopTools_IndexedMapOfShape& innerouterFaces,
		const std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair);

	TopoDS_Compound MakeGroups(const Face2CoedgesMap& faceCoedgesMap,
		const TopTools_IndexedMapOfShape& innerouterFaces,
		const std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair);

	void MatchSeamEdges(const std::vector<ContinueEdgePtr>& unmatched, std::vector<ContinueEdgePtr>& unmatchedAfterSeamEdgeMatch, CoedgeInfo& CCMatchedMap, double tol);

	// 辅助函数：合并face中的重复顶点
	TopoDS_Face MergeSameVerticesInFace(const TopoDS_Face& face);

	// 原有的face处理函数
	TopoDS_Face MergeSameVerticesForFace(const TopoDS_Face& face);

	TopoDS_Compound MergeSameVerticesForFaces(const TopoDS_Shape& shape);

	bool bPointsEqual(const gp_Pnt& p1, const gp_Pnt& p2, const double tol);

	void GetSoftVertices(const TopoDS_Shape& shape, double tol, TopTools_IndexedMapOfShape& softVertices, EdgeSplitMap& edgeSplitMap);

	void SplitEdgeByParams(const TopoDS_Edge& edge, const std::vector<double>& params, std::vector<TopoDS_Edge>& splited);

	TopoDS_Face SplitFaceEdges(const TopoDS_Face& face, const EdgeSplitMap& splitParams);

	TopoDS_Shape SplitShapeEdges(const TopoDS_Shape& shape, const EdgeSplitMap& splitParams);

	ContinueEdgePtr GetTwinContinueEdge(const ContinueEdgePtr& continueEdgePtr, const CoedgeInfo& edges);

	TopoDS_Shape GetBorderFacesForExtrudeShape(const TopoDS_Shape& shape);

	void FixEdgeInfo(const CoedgeInfo& coedgeInfoVector, const TopTools_IndexedMapOfShape& reversedFaces,
		const TopTools_IndexedDataMapOfShapeShape& faceMapping,
		CoedgeInfo& coedgeInfoVectorAfterFaceOritationAdjust);

	void GetFace2ContinueEdges(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfoVector, Face2CoedgesMap& face2Coedges, TopTools_IndexedMapOfShape& allFaces);

	TopTools_IndexedMapOfShape GetReversedFaces(const TopoDS_Shape& shape, const TopTools_IndexedMapOfShape& faces);

	bool FindCommonDirection(const TopoDS_Shape& shape, gp_Vec& dir,
		const Standard_Real theAngularTolerance = Precision::Angular());

	/**
	 * @brief 判断点集是否共面
	 * @param points 离散点集
	 * @param normal 如果共面, 返回平面法向
	 * @param squareTol 距离容差
	 * @return -1 : 非共面, 0: 点集退化, 1: 共面
	 */
	int IsPlanerCurve(const std::vector<gp_Pnt>& points, gp_Vec& normal, double squareTol = 1.0e-4);

	/**
	 * @brief 
	 * @param curve 待采样曲线 
	 * @param uMin 曲线参数域
	 * @param uMax 曲线参数域
	 * @param normal 返回共面法向
	 * @param squareTol 距离容差
	 * @return -1 : 非共面, 0: 点集退化, 1: 共面
	 */
	int IsPlanerCurve(const Handle(Geom_Curve)& curve, double uMin, double uMax, gp_Vec& normal,
		double squareTol = 1e-4);

	bool SampleContinue(const ContinueEdgePtr& continueEdge, int sampleNum, std::vector<gp_Pnt>& samplePnts);

	bool HasVRotationalAxis(const Handle(Geom_Surface)& surface, double uMin, double uMax, double vMin, double vMax, gp_Vec& direction);

	bool HasURotationalAxis(const Handle(Geom_Surface)& surface, double uMin, double uMax, double vMin, double vMax, gp_Vec& direction);

	bool HasSphereNearExtrudeDir(const TopoDS_Face& face, const gp_Vec& sphereNormal, gp_Vec& dir, double tol, double squareTol);

	bool GetSpecificSurfExtrudeDir(const Handle(Geom_Surface)& surf, gp_Vec& dir);

	TopoDS_Shape RemoveInnerLoopFaces(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfos, double thin);

	bool GetFaceInnerPt(const TopoDS_Face& face, double& u, double& v);

	void GetConnectFaceSeamEdges(const std::vector<std::pair<ContinueEdgePtr, bool>>& coedges, const CoedgeInfo& coedgeInfos, std::vector<ContinueEdgePtr>& seamEdges);

	void GetConnectFaceSeamEdges(const TopoDS_Face& face, const CoedgeInfo& coedgeInfos, std::vector<ContinueEdgePtr>& seamEdges);

	int SolidInOutStatus(const TopoDS_Solid& solid, bool isSingle, const gp_Vec& dir);

	gp_Cylinder GetAproxCylinder(const TopoDS_Shape& shape, const gp_Vec& dir, double& xmin, double& ymin, double& zmin, double& xmax, double& ymax, double& zmax);

	int ShapeInOutStatusByMesh(const TopoDS_Shape& shape, const gp_Vec& dir);

	int ShapeInOutStatusByMesh(const TopoDS_Shape& shape, const gp_Vec& dir, const gp_Cylinder& cylinder, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax);

private:  // un used api

	void GetEdgesSplitInterval(const TopoDS_Shape& shape, EdgeSplitMap& edgeSplitMap);

	Standard_Real DistancePointToEdge(const gp_Pnt& point, const TopoDS_Edge& edge, const Standard_Integer nbPnts = 10);

	// 计算点到Face边界的最小距离
	Standard_Real DistancePointToFaceBoundary(const gp_Pnt& point, const TopoDS_Face& face, double thinDist);

	double FindMaxDistancePointOnFace(const TopoDS_Face& face, double thinDist, double deflection = 0.01);

	bool IsConnectEdgeSmooth(const TopoDS_Shape& shape1, const TopoDS_Shape& shape2, const TopoDS_Vertex& v);

	bool GetMutualNearestPt(const std::shared_ptr<KDTree>& tree, const point_t& ptA, const pointVec& ptsVec, const double tol, point_t& ptNearest);

	size_t GetVertexCount(const std::vector<std::pair<gp_Pnt, int>>& vertexCount, const gp_Pnt& pt, double tol);

	size_t GetVertexCount(const std::vector<std::pair<gp_Pnt, std::vector<TopoDS_Edge>>>& vertexCount, const gp_Pnt& pt);

	bool CanMerge(const ContinueEdgePtr& continue1, const ContinueEdgePtr& continue2, bool checkEnd);

	const ContinueEdgePtr MergeEdges(const std::vector<ContinueEdgePtr>& continuesToMerge);

	std::vector<ContinueEdgePtr> MergeGroup(const std::vector<ContinueEdgePtr>& group);

	const std::vector<ContinueEdgePtr> MergeContinuousEdges(const std::vector<ContinueEdgePtr>& continues);

	bool MatchMergedEdges(const std::shared_ptr<KDTree>& tree, const std::vector<ContinueEdgePtr>& FirstUnMatched,
		CoedgeInfo& CCMatchedMap, std::vector<ContinueEdgePtr>& LastUnMatched, double tol, TopoDS_Compound& unMatchedEdges);

	double CDistance(const TopoDS_Face& face);

	double GetMaxLength(const TopoDS_Face& face, const std::vector<std::pair<gp_Pnt, int>>& vertexCount);

	double GetMaxLength(const TopoDS_Face& face, const std::vector<std::pair<gp_Pnt, std::vector<TopoDS_Edge>>>& vertexCount);

	double GetPerimeter(const TopoDS_Shape& comp);

	TopoDS_Face FixFaceWiresOrder(const TopoDS_Face& face);

	TopoDS_Shape RemoveShortEdges(const TopoDS_Shape& originalShape, Standard_Real minLength);

	double CalculateAreaExteriorOnly(const TopoDS_Face& face);

	void ComputeFaceNormals(
		const TopoDS_Shape& shape,
		gp_Dir& normalStart,    // at (umin,vmin)
		gp_Dir& normalMiddle,   // at (u50,v50)
		gp_Dir& normalEnd);

	double GetMinEdgeLength(const CoedgeInfo& coedgeInfos);

	TopoDS_Compound CreateCompoundFromFacePairs(
		const std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>>& facePairs);

	TopoDS_Shape MakeGroups(const CoedgeInfo& coedges, const TopTools_IndexedMapOfShape& inneroutFaces);

	TopoDS_Shape RemoveInnerLoopFaces(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfos, TopTools_IndexedMapOfShape& innerLoopFaces, double thin);

	bool IsFacesOnDifferentSide(const std::pair<ContinueEdgePtr, ContinueEdgePtr>& edge, const Face2CoedgesMap& map, double tol);

	bool JudgeTwoWileFaceBeCuttingFace(const TopoDS_Shape& face, const EdgeVectorPtr& faceCoedges, double thinDist);

	void GetContinuesEdges(const TopoDS_Shape& face, const CoedgeInfo& coedges, std::vector<std::pair<ContinueEdgePtr, ContinueEdgePtr>>& faceCoedges);

	void AddAdjecentSmoothFaces(const TopoDS_Shape& face, const Face2CoedgesMap& faceCoedgesMap,
		TopTools_IndexedMapOfShape& faceGroup, double tol);

	TopoDS_Compound GetSmoothFaces(const TopoDS_Face& face, const Face2CoedgesMap& faceCoedgesMap, double tol);
	/**
	 * @brief 查找仅出现一次的连续区间（若有多个，仅返回第一个满足条件的区间）。
	 *
	 * @param arr      输入数组
	 * @param startIdx 通过引用返回的区间起始位置
	 * @param endIdx   通过引用返回的区间结束位置
	 * @return true    找到只出现一次的连续区间
	 * @return false   未找到只出现一次的连续区间
	 */
	bool findSingleRun(const std::vector<int>& arr, int& startIdx, int& endIdx);

	TopoDS_Shape GetHolesOtherSideFaces(const TopoDS_Shape& holes, const TopoDS_Shape& smoothFaces, const Face2CoedgesMap& faceCoedgesMap, double tol);

	bool IsExtrude(const TopTools_IndexedMapOfShape& faces, gp_Vec& dir);

	bool IsOuterFace(const TopoDS_Face& face, BRepClass3d_SolidExplorer& aSE);

	bool HasAo(const EdgeVectorPtr& edges);

	AdjecentType CheckFacesDistribution(const TopTools_IndexedMapOfShape& currentPatchAdjecentFaces,
		const TopoDS_Shell& outerBorderFaces,
		const TopoDS_Shell& innerBorderFaces);

	AdjecentType CheckFacesDistribution(const TopTools_IndexedMapOfShape& currentPatchAdjecentFaces,
		const TopTools_IndexedMapOfShape& outerBorderFaces,
		const TopTools_IndexedMapOfShape& innerBorderFaces);

	void GetOuterFaceByRayInterersect(const shape2shapeMap& map, const TopoDS_Face& benchmarkFace, const CoedgeInfo& coedgeInfos,
		TopTools_IndexedMapOfShape& visited, BRepClass3d_SolidExplorer& aSE);

	TopoDS_Shape GetVertex(const std::vector<std::pair<gp_Pnt, int>>& vertexCountA);

	void GetMayBumpFaces(const CoedgeInfo& coedges, TopTools_IndexedMapOfShape& ret, double thinDist);

	TopoDS_Shape CollectFuzzyFaces(const TopTools_IndexedMapOfShape& fuzzyFaces,
		TopTools_IndexedMapOfShape& outerBorderFaces, TopTools_IndexedMapOfShape& innerBorderFaces,
		const Face2CoedgesMap& face2Coedges, double thinDist, double tol);

	void GetFace2ContinueEdges(const TopoDS_Shape& shape, const CoedgeInfo& coedgeInfoVector, Face2CoedgesMap& face2Coedges);

	void MakeTwoGroups(const TopoDS_Shape& patchs, TopTools_IndexedMapOfShape& outerBorderFaces, TopTools_IndexedMapOfShape& innerBorderFaces);

	std::vector<std::vector<ContinueEdgePtr>> PopLoop(const std::vector<int>& loop, const std::vector<ContinueEdgePtr>& totalContinues, double squareTol);

	void DiscretizeCurve(const Handle(Geom_Curve)& curve, double start, double end, int numPoints, std::vector<PointWithParam>& rets);

	void GetEdgePointsAndNormals(const TopoDS_Edge& edge,
		const TopoDS_Face& face,
		const int discreteNum,
		std::vector<gp_Pnt>& pnts, std::vector<gp_Dir>& vecs);

	void GetPointsAndNormalFromWires(std::vector<std::vector<ContinueEdgePtr>>& wires,
		const int discreteNum,
		std::vector<std::vector<gp_Pnt>>& pnts, std::vector<std::vector<gp_Dir>>& normals);

	bool IsCircularCurve(const std::vector<gp_Pnt>& points, gp_Pnt& center, gp_Vec& normal, double squareTol = 1.0e-4);

	bool IsCircularCurve(const Handle(Geom_Curve)& curve, double uMin, double uMax, gp_Pnt& center,
		double squareTol = 1e-4);

	// 辅助函数：解 (A - lambda * I) * x = 0，并把结果存到 vecOut
	// 若出现严重数值问题(如pivot很小)，返回false
	bool solveEigenVector3x3(const std::array<double, 9>& A,
		double lambda,
		std::array<double, 3>& vecOut);

	// 3x3特征值分解：
	// - 若成功，返回true，并在 outEigenPairs 中填入 3 个条目(对复根仅填实部虚部)
	// - 若出现数值问题(如奇异，NaN等)，返回false
	bool EigenDecompose3x3(const std::array<double, 9>& A,
		std::vector<EigenPair>& outEigenPairs);

	bool GetFaceExtrudeDirByPCA(const TopoDS_Face& face, gp_Vec& dir);

	bool GetExtrudeDirByPtsPCA(const TopoDS_Shape& shape, gp_Vec& dir);

	void GetRelatePatch(const TopoDS_Face& face, const std::vector<std::pair<TopoDS_Shell, TopoDS_Shell>>& mutuallyExclusivePair, TopTools_IndexedMapOfShape& curPatches, TopTools_IndexedMapOfShape& otherPatches);

	bool IsShellConnect(const TopoDS_Shell& shell1, const TopoDS_Shell& shell2, const Face2CoedgesMap& face2Coedges);

	void ClassifyFuzzyFacesByInOutGroups(TopoDS_Shape& outerBoundFaces, TopoDS_Shape& innerBoundFaces, const TopoDS_Compound& fuzzyUnRecognized,
		const Face2CoedgesMap& face2Coedges);
};

#endif