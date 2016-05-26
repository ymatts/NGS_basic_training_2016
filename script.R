#�K�v�ȃp�b�P�[�W�̓ǂݍ���
library(systemPipeR)
library(systemPipeRdata)

#�ψكR�[���p�̃e���v���[�g���[�N�t���[���쐬
#genWorkenvir(workflow="varseq")

#R�̍�ƃf�B���N�g����varseq�ɕύX
setwd("varseq")

#�K�v�Ȋ֐��̓ǂݍ���
source("systemPipeVARseq_Fct.R")

#targets�t�@�C���̃p�X���擾
targetspath <- "./targetsPE.txt"
#targets�t�@�C����ǂݍ���
targets <- read.delim(targetspath, comment.char = "#")[,1:5]
targets

#param�t�@�C���̃p�X���擾
parampath <- "./param/bwa.param"
#param�t�@�C����ǂݍ���
param <- read.delim(parampath, comment.char = "#")
param

#
args <- systemArgs(sysma=parampath, mytargets=targetspath)
names(args)

targetsin(args)
targetsout(args)
targetsheader(args)
modules(args)

#�K�v�Ȋ֐��̓ǂݍ���
source("systemPipeVARseq_Fct.R")
#args�̐ݒ�
args <- systemArgs(sysma="param/trimPE.param", mytargets="targetsPE.txt")[1:4]
#�e���[�h�z��ɂ����āAQUALITY�l<=20�𖞂���������v�Z���J�E���g�B���̃J�E���g����0�ȉ��̂��̂�����I��
#���[�h�z��Ɋ܂܂��QUALITY�l�����ׂ�20�ȏ�̔z�񂾂���I��
preprocessReads(args=args, Fct="filterFct(fq, cutoff=20, Nexceptions=0)", batchsize=100000)
targets_PEtrim.txt�ɏ�������
writeTargetsout(x=args, file="targets_PEtrim.txt", overwrite=TRUE)

#args�̐ݒ�
args <- systemArgs(sysma="param/bwa.param", mytargets="targets_PEtrim.txt")
#fastq�t�@�C������N�I���e�B���v���Z�o
fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
#result�f�B���N�g���̉���fastqReport.pdf�t�@�C���Ɍ��ʂ��o��
pdf("./results/fastqReport.pdf", height=18, width=4*length(fqlist))
seeFastqPlot(fqlist)
dev.off()

#args�̐ݒ�
args <- systemArgs(sysma="param/bwa.param", mytargets="targets_PEtrim.txt")
#
read_statsDF <- alignStats(args=args)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

#�K�v�p�b�P�[�W�̃C���X�g�[��
library(gmapR)
#gsnap�p�̃C���f�b�N�X�t�@�C�����쐬
gmapGenome <- GmapGenome("/sshare2/home/lect-1/BT2016/varseq/data/tair10.fasta",
directory="data",name="gmap_tair10chr",create=TRUE)
#args�̐ݒ�
args <- systemArgs(sysma="param/gsnap.param", mytargets="targets_PEtrim.txt")
#gsnap�̃p�����[�^�̐ݒ�
p <- GsnapParam(genome=gmapGenome, unique_only=TRUE, molecule="DNA", max_mismatches=3)
#gsnap�ɂ��A���C�����g
for(i in seq(along=args)) {
  o <- gsnap(input_a=infile1(args)[i], input_b=infile2(args)[i], params=p, output=outfile1(args)[i])
  print(c(infile1(args)[i],infile2(args)[i]))
}
file.exists(outpaths(args))
writeTargetsout(x=args, file="targets_gsnap_bam.txt", overwrite=TRUE)

#�K�v�ȃp�b�P�[�W�̓ǂݍ���
library(gmapR)
library(VariantTools)
#args�̐ݒ�
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

#�K�v�ȃp�b�P�[�W�̓ǂݍ���
library(VariantAnnotation)
library(BBmisc)
args <- systemArgs(sysma="param/filter_vartools.param", mytargets="targets_vartools.txt")
filter <- "(values(vr)$n.read.pos.ref + values(vr)$n.read.pos) >= 2 & (values(vr)$n.read.pos / (values(vr)$n.read.pos.ref + values(vr)$n.read.pos) >= 0.8)"
filterVars(args, filter, varcaller="vartools", organism="A. thaliana")
writeTargetsout(x=args, file="targets_vartools_filtered.txt", overwrite=TRUE)

#
#�K�v�ȃp�b�P�[�W�̓ǂݍ���
library("GenomicFeatures")
args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
#�A�m�e�[�V�����ɗp����f�[�^�x�[�X���w��
txdb <- loadDb("./data/tair10.sqlite")
#���t�@�����X�z��t�@�C�����w��
fa <- FaFile("/sshare2/home/lect-1/BT2016/varseq/data/tair10.fasta")
# �A�m�e�[�V�����̎��s
variantReport(args=args, txdb=txdb, fa=fa, organism="A. thaliana")

# �e�T���v���� nonsynonymous �ψق̕��������o������
combineDF <- combineVarReports(args=args, filtercol=c(Consequence="nonsynonymous"))
# ���������t�@�C����combineDF_nonsyn_vartools.xls�ɏ����o��
write.table(combineDF, "./results/combineDF_nonsyn_vartools.xls", quote=FALSE, row.names=FALSE, sep="\t")

args <- systemArgs(sysma="param/annotate_vars.param", mytargets="targets_vartools_filtered.txt")
# �ψكR�[���̃T�}���[��variantStats_vartools.xls�ɏ����o��
write.table(varSummary(args), "./results/variantStats_vartools.xls", quote=FALSE, col.names = NA, sep="\t")
