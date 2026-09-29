# oligopool-primer-design

[English](README.md) | **简体中文**

为寡核苷酸池（oligopool）合成设计**安全的通用引物结合位点**——适用于 MPRA、CRISPR 文库、barcode 池等任何需要用一对引物扩增整个池的混合寡核苷酸订单。

配套脚本 `build_oligos.py` 把选定的侧翼序列接到 payload 两端、生成可直接送合成的完整序列，并在写出前对**即将提交的 pool** 重新校验引物安全性。

## 问题背景

混合寡核苷酸订单（Twist、IDT 等）通常按如下结构合成：

```
[引物F结合序列] + payload（可变区） + [引物R结合序列（反向互补）]
```

这样整个池用一对 PCR 引物即可扩增。两条侧翼序列不能在 pool 内任何位置发生错配起扩——哪怕只有一个 3′ 端匹配到某条 payload，都会造成脱靶扩增、扭曲文库比例。

`oligopool-primer-design.py` 随机采样候选引物，只保留对 **pool 全部序列（含两条链）** 都安全的那些。

## 安全标准

候选引物必须同时满足：

1. GC 含量在 `[gc-min, gc-max]` 内（默认 40–60%）
2. 无同碱基连续 ≥ `max-homopolymer`（默认 5，如 AAAAA）
3. 无纯 G/C 连续 ≥ `max-gc-run`（默认 5）
4. 3′ 端 tail-k-mer **不出现**在 pool 任何序列（含反向互补链）中——防错配起扩的最关键一条
5. 全长任意 inner-k-mer 不出现在 pool 中
6. 任意两条已接受引物之间、以及引物自身的 3′ 端无互补（`dimer-k`，默认 8，设 0 关闭）——防止挑一对引物时出现 primer-dimer

pool 的反向互补链始终自动加入，因此这里通过的引物在两个方向上都安全。

> 注：经典做法"分别检查 11-mer 和 12-mer"是冗余的——任何 12-mer 命中必然伴随其子 11-mer 命中。本工具采用单一内部 k-mer 规则（自动选择，见下），等效于两者中更严格的那个。

## k 阈值如何确定

tail-k / inner-k **完全自动、不可手动指定**——程序总是为你的 pool 选择可行的最小组合：从最严格的 10 / 11 出发，用 15 万条随机候选做预试验、对实际 pool 实测通过率，不行就 +1 逐步放宽（上限为 min(21, 引物长度)），直到有候选通过。每次运行都会打印选定的组合。这保证每个设计都用其 pool 所允许的最严格标准，无需猜测任何参数。

参考值（双链总碱基数；来自随机序列与基因组来源测试 pool 的实测）：

| 双链总碱基数 | 典型可行最小 tail-k / inner-k |
|---|---|
| 2,000,000 bp 以下 | 10 / 11 |
| 2,000,000 – 16,000,000 bp | 10 / 11 – 11 / 12 |
| 16,000,000 – 60,000,000 bp | 11 / 12 – 13 / 14 |
| 60,000,000 bp 以上 | 13 / 14 或更高 |

基因组来源的 pool 通常比同样大小的随机序列 pool 允许更小的 k：其 k-mer 分布受 GC 偏向影响，而候选引物 GC 居中，恰好落在稀有 k-mer 区域。若到上限仍无可行组合，程序会明确报告而不是死循环。

## 环境要求

- Python 3.6+（仅标准库）
- `numpy` **可选**——大 pool 时 k-mer 集合构建显著加速（未安装自动回退纯 Python）

## 快速上手

```bash
# 默认标准生成 20 条候选引物
# (k 阈值自动选择: 取当前 pool 可行的最小组合)
python3 oligopool-primer-design.py --fasta my_pool.fa

# 从中挑一对, 比如 P03 作正向、P11 作反向
#   5′ 侧翼 = P03
#   3′ 侧翼 = P11        (其反向互补才是实际的反向引物)

# 可复现运行(单线程/多线程结果完全一致)
python3 oligopool-primer-design.py --fasta my_pool.fa --seed 42

# pool 更新后, 重新校验之前设计的引物是否仍然安全
python3 oligopool-primer-design.py --fasta my_pool_v2.fa --check safe_primers.txt
```

输出：`safe_primers.txt` 每行一条引物（方便管道处理），另附带 GC 与 Wallace Tm 的 `.tsv` 明细。`--check` 在任一引物不通过时以非零码退出，可直接接入 CI/Snakemake。

## 生成送合成的完整序列

`build_oligos.py` 接收 payload FASTA 和选定的两条引物，输出可直接下单的完整序列：

```bash
python3 build_oligos.py --fasta my_pool.fa \
    --primer-f CTCGATGATGAAAACCGTCT --primer-r TCCCGACTAAGCCCATGGAT
```

产出（前缀 `--out`，默认 `oligo_for_synthesis`）：

- `out.fa` — `[引物F] + payload + [引物R位点]`，保留原始名称，80 列换行
- `out.csv` — `name,sequence` 两列，合成厂商订单格式
- `out.pcr_primers.txt` — 实际使用的 PCR 引物对（F 原样，R 为 3′ 侧翼的反向互补），省得拿到文件的人自己换算方向

写出前会**对即将提交的 pool 重新校验两条引物**（3′ 端 k-mer、内部 k-mer、二聚体；k 同样自动选择）。若引物已不安全——比如 pool 在引物设计之后扩充过——则拒绝写出并以非零码退出（`--force` 强制写出，`--no-verify` 跳过校验）。同时会剔除含非 ACGT 字符的序列（无法合成），并在最终 oligo 超过 `--max-len`（默认 300 nt）时告警。

方向约定：

- `--primer-f` 加在 5′ 端，本身就是 PCR 正向引物。
- `--primer-r` 按原样加在 3′ 端；PCR 反向引物是它的反向互补。如果你手上拿的已经是反向引物本身，加 `--primer-r-is-rc`，程序会自动取反向互补后再拼接。

## 参数 — `oligopool-primer-design.py`

| 参数 | 默认值 | 含义 |
|---|---|---|
| `--fasta` | （必填） | oligopool 序列，FASTA 格式 |
| `--out` | fasta 同目录 `safe_primers.txt` | 输出文件 |
| `-n, --n-primers` | 20 | 生成引物条数 |
| `--threads` | 20 | 并行进程数（1 = 单线程） |
| `--primer-len` | 20 | 引物长度 |
| `--gc-min / --gc-max` | 0.40 / 0.60 | GC 含量窗口 |
| `--max-homopolymer` | 5 | 同碱基连续达到该值即拒绝 |
| `--max-gc-run` | 5 | 纯 G/C 连续达到该值即拒绝 |
| `--dimer-k` | 8 | 引物间 3′ 端互补检查长度（0 = 关闭） |
| `--gc-clamp` | 0 | 要求 3′ 端 5 nt 内至少 N 个 G/C（0 = 关闭） |
| `--seed` | 无 | 随机种子，完全可复现 |
| `--max-attempts` | 2,000,000 | 每条引物最大尝试次数（防死循环） |
| `--check FILE` | — | 校验已有引物对当前 pool 的安全性 |

## 参数 — `build_oligos.py`

| 参数 | 默认值 | 含义 |
|---|---|---|
| `--fasta` | （必填） | payload 序列，FASTA 格式 |
| `--primer-f` / `--primer-r` | （必填） | 两侧引物结合序列 |
| `--primer-r-is-rc` | 关 | 给定的 `--primer-r` 已是反向引物本身，先取反向互补再拼接 |
| `--out` | fasta 同目录 `oligo_for_synthesis` | 输出前缀（`.fa` / `.csv` / `.pcr_primers.txt`） |
| `--no-verify` | 关 | 跳过引物-pool 安全性校验 |
| `--force` | 关 | 校验不通过也强制写出 |
| `--dimer-k` | 8 | 校验用二聚体检查长度（0 = 关闭） |
| `--max-len` | 300 | 最终 oligo 超过该长度时告警 |

## 性能

42,000 条 160 bp oligo（6.7 Mb pool，含反向互补 13.4 Mb）：常规实验室服务器 + numpy 下**总计约 10 秒**——k-mer 集合用向量化滚动 2-bit 编码构建、经 fork 写时复制共享给各 worker，k 自动选择含 15 万条预试验验证，引物由并行拒绝采样生成。

## 许可证

MIT——见 [LICENSE](LICENSE)。
