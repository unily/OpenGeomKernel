
#ifndef _CUTTINGRAILSOPTION_HeaderFile
#define _CUTTINGRAILSOPTION_HeaderFile

#include "../COMMON/Extension_Export.hxx"
#include <cstdint>

/**
 * @brief SingleBorderType：定义模型的单壁类型
 *
 *  - None：不确定，可能是单壁也可能是双壁
 *          （算法内部推断：当模型壁厚 > 30 时被视为单外壁）
 *  - Outer：仅外壁
 *  - Inner：仅内壁
 */
EXTENSION_EXPORT enum class SingleBorderType : uint8_t
{
    None = 0,  //!< 不确定，可能是单壁也可能是双壁
    Outer = 1,  //!< 仅外壁
    Inner = 2   //!< 仅内壁
};

/**
 * @brief CuttingRailsOption：配置刀轨识别选项
 */
EXTENSION_EXPORT struct CuttingRailsOption
{
    /**
     * @brief 指示输入模型是否为等径模型, 不确定是等径还是非等径时填false
     * 默认为 false。
     */
    bool bExtrudeShape{ false };

    /**
     * @brief 是否需要切割面，5轴切割为true，4轴切割为false
     * 默认为 false。
     */
    bool bNeedCuttingFaces{ false };

    /**
     * @brief 输入模型的单壁类型
     * 默认为 SingleBorderType::None。
     */
    SingleBorderType singleBorderType{ SingleBorderType::None };

    /**
     * @brief 默认构造函数
     */
    CuttingRailsOption() = default;

    /**
     * @brief 使用自定义值初始化所有成员
     * @param extrude    是否为等径模型
     * @param borderType 模型的单壁类型
     */
    CuttingRailsOption(bool needCuttingFaces, bool extrude, SingleBorderType borderType)
        : bNeedCuttingFaces(needCuttingFaces), bExtrudeShape(extrude), singleBorderType(borderType)
    {
    }
};

#endif