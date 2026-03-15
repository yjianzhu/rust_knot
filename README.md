# rust_knot

Rust 实现的纽结拓扑分析库 — 通过 Alexander 多项式识别纽结类型并定位最小纽结核心。

本项目是对 C++ 版本（`knottype.cpp`, `knotsize.cpp`, `hull.cpp`, `my_function.cpp`）的完整重写，消除了 GiNaC 符号数学依赖，修复了已知 bug，并提供线程安全的纯函数式设计。

## 算法流程

```
输入点序列 (3D 坐标)
    │
    ▼
KMT 简化 (去除不影响拓扑的冗余点)
    │
    ▼
凸包端点延伸 (开链: 将首尾推离纽结区域)
    │
    ▼
O(n²) 交叉点检测 (XY 投影 + Z 序比较)
    │
    ▼
Alexander 矩阵构造 (Z[t] 系数)
    │
    ▼
Bareiss 无分数行列式 (替代 GiNaC)
    │
    ▼
多项式查表匹配
    │
    ▼
纽结类型 ("3_1", "4_1", ...)
    │
    ▼
二分搜索最小纽结核心 (knotsize)
    │
    ▼
KnotCoreResult { left, right, size }
```

## 项目结构

```
rust_knot/
├── Cargo.toml
├── src/
│   ├── lib.rs               # crate 入口，re-export 公开 API
│   ├── main.rs              # CLI 入口（直接参数模式）
│   ├── config.rs             # KnotConfig — 统一配置结构体
│   ├── batch.rs              # 批处理: process_frame / process_frames_parallel / process_frames_streaming
│   ├── point.rs              # Point3 = [f64; 3], EPSILON = 1e-7
│   ├── error.rs              # KnotError 枚举 (Io, PolynomialParse, DataParse, HullFailed, NotFound, EmptyChain)
│   ├── polynomial.rs         # Polynomial<i64> 算术 + Bareiss 行列式 (~600 行，替代 GiNaC)
│   ├── geometry.rs           # 交叉检测、三角形相交测试、法线计算、坐标轴重定向
│   ├── alexander_table.rs    # Alexander 多项式查找表 (内嵌 ≤9 crossings + 外部文件合并)
│   ├── kmt.rs                # KMT 简化算法 (开链 + 环链)
│   ├── hull.rs               # 凸包端点延伸 (基于 chull crate)
│   ├── knottype.rs           # 纽结类型识别核心: 交叉点 → Alexander 矩阵 → 多项式 → 查表
│   ├── knotsize.rs           # 二分搜索纽结核心定位
│   └── io.rs                 # XYZ / LAMMPS 格式读写 + XyzFrameIter 惰性帧迭代器
├── tests/
│   └── integration.rs        # 6 个端到端集成测试
└── .github/
    └── workflows/
        ├── ci.yml            # GitHub Actions: push/PR 代码检查
        └── release.yml       # GitHub Actions: tag 触发三平台发布
```

## 依赖

| crate | 版本 | 用途 |
|-------|------|------|
| `chull` | 0.2 | 3D 凸包 (QuickHull 算法) |
| `thiserror` | 2 | 错误类型派生 |
| `rayon` | 1.10 | 帧间并行批处理 |
| `approx` | 0.5 | 测试用浮点比较 (dev-dependency) |

全部为纯 Rust 依赖，无 C/FFI 绑定，支持零配置交叉编译。

**关键设计决策**: 不使用任何符号数学库。Alexander 矩阵的元素只有 `{0, 1, -1, t, 1-t}`，因此我们实现了轻量的 `Polynomial { coeffs: Vec<i64> }`（约 600 行），配合 Bareiss 无分数消元行列式算法，完全在 `Z[t]` 整数多项式环内完成计算。

## 配置 (`KnotConfig`)

所有超参数和模式标志统一管理在 `KnotConfig` 结构体中:

| 字段 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `is_ring` | `bool` | `false` | 开链 (`false`) 或环链 (`true`) |
| `faster` | `bool` | `false` | 启用 KMT 简化（加速，不影响结果） |
| `debug` | `bool` | `false` | 输出调试信息到 stderr |
| `hull_plane_epsilon` | `f64` | `5e-3` | 凸包面判定阈值（越大容忍度越高） |
| `extend_factor` | `f64` | `100.0` | 端点外延缩放因子（越大端点推得越远） |
| `num_rotations` | `u32` | `4` | 环链模式旋转搜索次数（越多越精确但越慢） |

`EPSILON = 1e-7`（`point.rs`）为全局几何零判定阈值，用于交叉检测和法线计算，作为底层数值常量不纳入 `KnotConfig`。

## 编译和测试

```bash
cd rust_knot

# 编译
cargo build --release

# 运行全部 51 个测试 (45 单元 + 6 集成)
cargo test

# 编译 CLI 工具
cargo build --release
```

## 使用

### CLI 工具

```bash
# 最简用法: 内嵌表 (≤9 crossings)，无需外部文件
cargo run -- input.xyz
# 或 release 二进制
./target/release/rust_knot input.xyz

# 环链模式
./target/release/rust_knot input.xyz --ring

# 追加外部表 (>9 crossings)
./target/release/rust_knot input.xyz --table extended_table.txt

# 多帧文件 + 指定批大小 + 自定义输出路径
./target/release/rust_knot trajectory.xyz --ring --batch 128 --threads 8 --output result.txt

# 全部参数
./target/release/rust_knot <xyz_file> [target_type] [--table <path>] [--ring] [--fast] [--debug] [--output <path>] [--batch <size>] [--threads <n>]
```

输出文件 `knot_index.txt`:
```
# frame	knottype	knot_start	knot_end	knot_size
0	3_1	114	145	32
1	3_1	107	158	52
2	3_1	108	171	64
...
```

### 作为库

```rust
use rust_knot::{AlexanderTable, KnotConfig, get_knottype, find_knot_core};

// 方式一: 使用内嵌表 (≤9 crossings，无需外部文件)
let table = AlexanderTable::builtin();

// 方式二: 内嵌表 + 外部文件合并 (覆盖 ≥10 crossings)
let table = AlexanderTable::builtin_with_file("extended_table.txt")?;

// 方式三: 仅使用外部文件
let table = AlexanderTable::from_file("full_table.txt")?;

// 配置
let config = KnotConfig {
    faster: true,
    is_ring: false,
    ..KnotConfig::default()
};

// 识别纽结类型
let knot_type = get_knottype(&points, &table, &config)?;  // e.g. "3_1"

// 定位最小纽结核心
let core = find_knot_core(&points, &knot_type, &table, &config)?;
println!("核心区间: [{}, {}], 大小: {}", core.left, core.right, core.size);
```

#### 多帧批处理

```rust
use rust_knot::batch::{process_frames_parallel, process_frames_streaming};
use std::io::BufReader;
use std::fs::File;

// 方式一: 全部加载 + 并行处理 (小文件)
let results = process_frames_parallel(&frames, &table, &config, None);

// 方式二: 流式批处理 (大文件，内存受控)
let reader = BufReader::new(File::open("huge_trajectory.xyz")?);
let total = process_frames_streaming(
    reader, &table, &config, None,
    Some(128),           // batch_size
    |batch_results| {    // 每批回调
        for r in batch_results {
            println!("frame {}: {}", r.frame, r.knot_type);
        }
    },
)?;
```

## Alexander 多项式表

### 内嵌表

程序内嵌了 crossing number ≤ 9 的全部 86 条 Alexander 多项式（unknot + 3₁ 到 9₄₉），覆盖绝大多数常见纽结类型。使用 `AlexanderTable::builtin()` 即可，无需任何外部文件。会使用正负两份存入哈希表，确保正反转都能识别。

### 外部扩展表

如需识别 ≥ 10 crossings 的纽结，通过 `--table` 参数或 `builtin_with_file()` 加载外部表。格式:

```
knot_name	polynomial
10_1	4-9*t+4*t^2
10_2	-2+7*t-2*t^2
...
```

内嵌表与外部表自动合并去重，同一多项式对应多种纽结时按 crossing number 排序（最简者优先）。

## 跨平台发布

GitHub Actions 自动构建三平台 release：

| 平台 | 产物 |
|------|------|
| Linux x86_64 | `rust_knot-linux-x86_64` |
| Windows x86_64 | `rust_knot-windows-x86_64.exe` |
| macOS ARM | `rust_knot-macos-arm64` |

发布新版本:
```bash
git tag v0.2.0
git push origin v0.2.0
# CI 自动编译 + 发布到 GitHub Releases
```

## 相对 C++ 版本的改进

### Bug 修复

| Bug | C++ 位置 | Rust 修复 |
|-----|----------|-----------|
| 字符比较判断纽结复杂度 | `knotsize.cpp:259,279` — `temp[0] > target[0]` | `parse_knot_name("10_1") → (10,1)` 数值比较 |
| 未知多项式返回空字符串 | `knottype.cpp:281` | 返回 `Err(KnotError::NotFound)` |
| NaN 输入导致 panic | `partial_cmp().unwrap()` | 使用 `total_cmp()` NaN 安全排序 |
| 全局 GiNaC mutex 阻塞并行 | `knottype.cpp:193` | 无共享可变状态，rayon 帧间并行 |
| 凸包失败静默忽略 | `hull.cpp catch(...)` | 返回 `Option`，debug 模式输出警告 |
| 查找表格式错误静默跳过 | 无校验 | strict 模式（默认）报错含行号，lenient 模式向后兼容 |

### 架构改进

- **无 GiNaC 依赖**: 整数多项式 + Bareiss 行列式完全替代符号计算库
- **内嵌多项式表**: ≤9 crossings 内建，零外部文件启动
- **流式批处理**: `XyzFrameIter` 惰性读取 + 分批 rayon 并行，内存占用恒定
- **线程安全**: 所有核心函数为纯函数，无全局状态
- **跨平台**: 纯 Rust 依赖，GitHub Actions 三平台自动发布
- **强类型错误处理**: `KnotError` 枚举覆盖所有错误路径

## 测试覆盖

51 个测试 (45 单元 + 6 集成):

| 模块 | 测试数 | 覆盖内容 |
|------|--------|----------|
| `polynomial` | 16 | 四则运算、行列式、解析、归一化、精确除法、相等/哈希 |
| `alexander_table` | 7 | 表解析、正负号查找、歧义处理、strict/lenient 模式、内嵌表、合并 |
| `geometry` | 4 | 交叉检测、法线计算、叉积 |
| `hull` | 3 | 端点延伸、退化几何、短链 |
| `kmt` | 2 | 直线简化、端点保留 |
| `knotsize` | 3 | 纽结名解析、复杂度比较、点序列旋转 |
| `io` | 4 | XYZ 读写、多帧读取、往返一致性 |
| `batch` | 2 | 单帧处理、并行结果有序性 |
| **集成测试** | **6** | trefoil ring、figure-eight open、unknot、knot core 定位、歧义表查询 |

## License

与上级项目保持一致。
