#include "CuttingRails.hxx"
#include "CuttingRailsImpl.hxx"
#include <iostream>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <ShapeUpgrade_RemoveLocations.hxx>

void CuttingRails::GetCuttingRails(const TopoDS_Shape& shape,
	TopoDS_Shape& outerBoundFaces,
	TopoDS_Shape& innerBoundFaces,
	TopoDS_Compound& rails,
	TopoDS_Compound& railsBoderFaces,
	TopoDS_Compound& twinRails,
	TopoDS_Compound& twinRailsCuttingFaces,
	TopoDS_Compound& unMatchedEdges,
	const CuttingRailsOption& option,
	const double tolerance)
{
	CuttingRailsImpl impl;

	bool hasface = impl.HasFace(shape);
	if (!hasface) {
		// std::cout << "has no face" << std::endl;
		return;
	}

	TopoDS_Shape comp = impl.copyShape(shape);
	ShapeUpgrade_RemoveLocations locationFix;
	locationFix.Remove(comp);
	comp = locationFix.GetResult();
	double tol = std::max(std::min(0.1, tolerance), 0.001);
	std::vector<std::vector<CuttingRailsImpl::ContinueEdgePtr>> loops;
	std::vector<std::vector<CuttingRailsImpl::ContinueEdgePtr>> twinLoops;
	impl.GetCuttingRails(comp, outerBoundFaces, innerBoundFaces, unMatchedEdges, loops, twinLoops, option, tol);

	// 输出loops
	BRep_Builder BB;
	if (rails.IsNull()) {
		BB.MakeCompound(rails);
	}

	if (railsBoderFaces.IsNull()) {
		BB.MakeCompound(railsBoderFaces);
	}

	if (twinRails.IsNull()) {
		BB.MakeCompound(twinRails);
	}

	if (twinRailsCuttingFaces.IsNull()) {
		BB.MakeCompound(twinRailsCuttingFaces);
	}

	for (auto& loop : loops)
	{
		TopoDS_Wire theWire;
		BB.MakeWire(theWire);
		for (auto& shp : loop)
		{
			for (const auto& e : shp->edges)
			{
				BB.Add(theWire, e);
				BB.Add(railsBoderFaces, shp->face);
			}
		}
		BB.Add(rails, theWire);
	}


	for (auto& loop : twinLoops)
	{
		TopoDS_Wire theWire;
		BB.MakeWire(theWire);
		for (auto& shp : loop)
		{
			if (!shp) {
				continue;
			}
			for (const auto& e : shp->edges)
			{
				BB.Add(theWire, e);
				BB.Add(twinRailsCuttingFaces, shp->face);
			}
		}
		BB.Add(twinRails, theWire);
	}

	return;
}

bool CuttingRails::GetExtrudeDirection(const TopoDS_Shape& shape, gp_Vec& dir) {
	CuttingRailsImpl impl;
	std::vector<double> faceAreas;
	TopoDS_Shape shapeAfterFaceSorted = impl.SortFaceByArea(shape, faceAreas);
	return impl.GetSortedShapeExtrudeDir(shapeAfterFaceSorted, faceAreas, dir);
}

bool CuttingRails::GetIsolatedEdges(const TopoDS_Shape& shape, TopoDS_Compound& isolatedEdges)
{
	BRep_Builder builder;
	if (isolatedEdges.IsNull()) {
		builder.MakeCompound(isolatedEdges);
	}

	TopExp_Explorer exp(shape, TopAbs_EDGE, TopAbs_FACE);
	for (; exp.More(); exp.Next()) {
		auto& edge = exp.Current();
		builder.Add(isolatedEdges, edge);
	}
	return true;
}
