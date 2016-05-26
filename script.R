#必要なパッケージの読み込み
library(systemPipeR)
library(systemPipeRdata)

#変異コール用のテンプレートワークフローを作成
#genWorkenvir(workflow="varseq")

#Rの作業ディレクトリをvarseqに変更
setwd("varseq")

#必要な関数の読み込み
source("systemPipeVARseq_Fct.R")

#targetsファイルのパスを取得
targetspath <- "./targetsPE.txt"
#targetsファイルを読み込み
targets <- read.delim(targetspath, comment.char = "#")[,1:5]
targets

#paramファイルのパスを取得
parampath <- "./param/bwa.param"
#paramファイルを読み込み
param <- read.delim(parampath, comment.char = "#")
param

#
args <- systemArgs(sysma=parampath, mytargets=targetspath)
names(args)

targetsin(args)
targetsout(args)
targetsheader(args)
modules(args)

#必要な関数の読み込み
source("systemPipeVARseq_Fct.R")
#argsの設定
args <- systemArgs(sysma="param/trimPE.param", mytargets="targetsPE.txt")[1:4]
#各リード配列において、QUALITY値<=20を満たす塩基数を計算しカウント。そのカウント数が0以下のものだけを選択
#リード配列に含まれるQUALITY値がすべて20以上の配列だけを選択
preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
targets_PEtrim.txtに書き込み
writeTargetsout(x=args, file="targets_PEtrim.txt", overwrite=TRUE)

#argsの設定
args <- systemArgs(sysma="param/bwa.param", mytargets="targets_PEtrim.txt")
#fastqファイルからクオリティ統計を算出
fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
#resultディレクトリの下のfastqReport.pdfファイルに結果を出力
pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
seeFastqPlot(fqlist)
dev.off()

#argsの設定
args <- systemArgs(sysma="param/bwa.param", mytargets="targets_PEtrim.txt")
#
read_statsDF <- alignStats(args=args)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

#必要パッケージのインストール
library(gmapR)
#gsnap用のインデックスファイルを作成
gmapGenome <- GmapGenome("/sshare2/home/lect-1/BT2016/varseq/data/tair10.fasta",
directory="data",name="gmap_tair10chr",create=TRUE)
#argsの設定
args <- systemArgs(sysma="param/gsnap.param", mytargets="targets_PEtrim.txt")
#gsnapのパラメータの設定
p <- GsnapParam(genome=gmapGenome, unique_only=TRUE, molecule="DNA", max_mismatches=3)
#gsnapによるアライメント
for(i in seq(along=args)) {
  o <- gsnap(input_a=infile1(args)[i], input_b=infile2(args)[i], params=p, output=outfile1(args)[i])
  print(c(infile1(args)[i],infile2(args)[i]))
}
file.exists(outpaths(args))
writeTargetsout(x=args, file="targets_gsnap_bam.txt", overwrite=TRUE)

#必要なパッケージの読み込み
library(gmapR)
library(VariantTools)
#argsの設定
args <- systemArgs(sysma="param/vartools.param", mytargets="targets_gsnap_bam.txt")
gmapGenome <- GmapGenome("/sshare2/home/lect-1/BT2016/varseq/data/tair10.fasta",
directory="data", name="gmap_tair10chr", create=FALSE)
tally.param <- TallyVariantsParam(gmapGenome, high_base_quality = 23L, indels = TRUE)
for(x in seq(along=args)){
  bfl <- BamFileList(infile1(args)[x], index=character())
  var <- callVariants(bfl[[1]], tally.param)
  sampleNames(var) <- names(bfl)
  writeVcf(asVCF(var), outfile1(args)[x], index = TRUE)
  print(outfile1(args)[x])
}
file.exists(outpaths(args))
writeTargetsout(x=args, file="targets_vartools.txt", overwrite=TRUE)

#必要なパッケージの読み込み
library(VariantAnnotation)
library(BBmisc)
args <- systemArgs(sysma="param/filter_vartools.param", mytargets="targets_vartools.txt")
filter <- "(values(vr)$n.read.pos.ref + values(vr)$n.read.pos) >= 2 & (values(vr)$n.read.pos / (values(vr)$n.read.pos.ref + values(vr)$n.read.pos) >= 0.8)"
filterVars(args, filter, varcaller="vartools", organism="A. thaliana")
writeTargetsout(x=args, file="targets_vartools_filtered.txt", overwrite=TRUE)

#
#必要なパッケージの読み込み
library("GenomicFeatures")
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
#アノテーションに用いるデータベースを指定
txdb <- loadDb("./data/tair10.sqlite")
#リファレンス配列ファイルを指定
fa <- FaFile("/sshare2/home/lect-1/BT2016/varseq/data/tair10.fasta")
# アノテーションの実行
variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana")

# 各サンプルの nonsynonymous 変異の部分を取り出し結合
combineDF <- combineVarReports(args=args, filtercol=c(Consequence="nonsynonymous"))
# 結合したファイルをcombineDF_nonsyn_vartools.xlsに書き出し
write.table(combineDF, "./results/combineDF_nonsyn_vartools.xls", quote=FALSE, row.names=FALSE, sep="\t")

args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
# 変異コールのサマリーをvariantStats_vartools.xlsに書き出し
write.table(varSummary(args), "./results/variantStats_vartools.xls", quote=FALSE, col.names = NA, sep="\t")

