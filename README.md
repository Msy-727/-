# STL 诊断与修复技能 · stl-diagnose-repair

针对 3D 打印切片软件（如 Raise3D ideaMaker）报「非法的边 / 非流形边」的 STL 模型，提供**拓扑+几何深度诊断 → 可视化报告 → 干净网格重建**的端到端工具链。

适用于 SolidWorks 建模 + ideaMaker 切片的工作流。

## 功能
- **深度诊断**：非流形边、共面重叠面、自相交、近重合顶点（裂缝隐患）、连通壳体数、欧拉特征（自动判定孔数/亏格）。
- **可视化报告**：单文件 HTML，平板件给截面轮廓预览、曲面件/壳体件给正交三视图 + 切片扫描，纯 SVG 无外部依赖。
- **网格重建**：从任意等厚拉伸件（支架/标牌/字母/垫片/齿条）自动提取轮廓并重建绝对干净的 STL，支持改厚度与多孔。
- **终端预览**：ASCII 三视图，无需图像库即可快速确认形状。

## 脚本（scripts/）
| 脚本 | 用途 |
|------|------|
| `deep_check.py <模型.stl>` | 6 项深度拓扑/几何检查 |
| `make_diag_report.py <输入.stl> <输出.html>` | 平板件截面诊断报告 |
| `make_view_report.py <输入.stl> <输出.html>` | 曲面件/壳体件三视图报告 |
| `ascii_view.py <输入.stl> [宽]` | 终端 ASCII 三视图 |
| `rebuild_extrude.py <输入.stl> <输出.stl> [厚度]` | 重建干净网格 / 改厚度 |

## 用法示例
```bash
python scripts/deep_check.py model.stl
python scripts/make_diag_report.py model.stl report.html
python scripts/rebuild_extrude.py model.stl clean.stl 5   # 重建为 5mm 厚
```

## 已知踩坑（详见 SKILL.md）
- 不能假设厚度沿 Z 轴（某些件沿 X 拉伸），需自动识别最短轴作为拉伸方向。
- 共面重叠不能用「顶点包含数」判定，否则把正常相邻面误判为非法；必须用多边形裁剪算真实交叠面积。
- 多孔桥接必须一次性构造单环，逐个插入会让外环边索引失效。
- 桌面中文名 STL 易被移动，`ls` 查不到，需用 Python 绝对路径定位。

## 依赖
仅 Python 3 标准库，无需 pip 安装任何包。
