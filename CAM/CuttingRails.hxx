
#ifndef _CUTTINGRAILS_HeaderFile
#define _CUTTINGRAILS_HeaderFile

#include "../COMMON/Extension_Export.hxx"
#include "./CuttingRailsOption.hxx"

class TopoDS_Shape;
class TopoDS_Compound;
class gp_Vec;

class CuttingRails {
public:
	/**
	 * @brief 提取模型的刀具轨迹
	 * @param [in] shape 待提取的模型，要求为薄壁模型
	 * @param [out] outerBoundFaces  识别出来的外壁面
	 * @param [out] innerBoundFaces  识别出来的内壁面
	 * @param [out] rails 识别出来的轨迹线, 由多个wire组成，在外壁面上的轨迹线
	 * @param [out] railsBoderFaces 轨迹线对应的壁面
	 * @param [out] twinRails 在切割面上的轨迹线，和rails是配对的。如果bNeedCuttingFaces是false, twinRails为空
	 * @param [out] twinRailsCuttingFaces 轨迹线对应的切割面。如果bNeedCuttingFaces是false, twinRailsCuttingFaces为空
	 * @param [in] option 刀轨识别的相关配置参数
	 * @param [in] tolerance 共边匹配容差
	 */
	EXTENSION_EXPORT static void GetCuttingRails(const TopoDS_Shape& shape,
		TopoDS_Shape& outerBoundFaces,
		TopoDS_Shape& innerBoundFaces,
		TopoDS_Compound& rails,
		TopoDS_Compound& railsBoderFaces,
		TopoDS_Compound& twinRails,
		TopoDS_Compound& twinRailsCuttingFaces,
		TopoDS_Compound& unMatchedEdges,
		const CuttingRailsOption& option,
		const double tolerance = 0.1);

	/**
	* @brief 提取拉伸方向
	* @param [in] shape 模型
	* @param [out] 拉伸方向
	* @param [out] isLicenseValid 当license在有效期, 返回true; 否则返回false
	*/
	EXTENSION_EXPORT static bool GetExtrudeDirection(const TopoDS_Shape& shape, gp_Vec& dir, bool& isLicenseValid);

	/**
	 * @brief 获取模型孤立边
	 * @param [in] shape 模型
	 * @param [out] isolatedEdges 孤立边 
	 * @param [out] isLicenseValid 当license在有效期, 返回true; 否则返回false
	 * @return 
	 */
	EXTENSION_EXPORT static bool GetIsolatedEdges(const TopoDS_Shape& shape, TopoDS_Compound& isolatedEdges, bool& isLicenseValid);
};

#endif