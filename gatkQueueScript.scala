import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.utils.commandline.{Output, Argument, Hidden}
import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.picard._
//import htsjdk.samtools.{SAMReadGroupRecord, SAMFileReader, SAMFileHeader}
import org.broadinstitute.gatk.utils.baq.BAQ.CalculationMode
import collection.JavaConversions._
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException
//import htsjdk.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.gatk.tools.walkers.indels.IndelRealigner.ConsensusDeterminationModel
//import org.broadinstitute.gatk.engine.arguments.ValidationExclusion
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.gatk.tools.walkers.genotyper.OutputMode
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypingOutputMode
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode
import org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode

class gatkpipeline extends QScript {
  qscript =>
  
  @Input(doc="bam file list. -- one per -I",
    fullName="input", shortName="I", required=true)
  var input: Seq[File] = _
  
  @Input(doc="The reference file for the bam files.", shortName="R")
  var reference: File = _
  
  @Input(doc="dbsnp",
    fullName="dbsnp", shortName="D", required=true)
  var dbSNP: Seq[File] = _ 
  
  @Argument(doc="Which genotyper you want to use: UnifiedGenotyper, or HaplotypeCaller, or UnifiedGenotyperAndHaplotypeCaller", 
    fullName="genotyper", shortName="genotyper", required=true)
  var caller: String = _
  
  @Argument(doc="split N Cigar Reads or not.",
    fullName="splitNCigarReads", shortName="splitNCigarReads", required=false)
  var splitNCigarReads = false
  
  @Argument(doc="filter reads with N cigar",
    fullName="filter_reads_with_N_cigar", shortName="filter_reads_with_N_cigar", required=false)
  var filter_reads_with_N_cigar = false


  @Argument(doc="do not use soft clipped bases.",
    fullName="dontUseSoftClippedBases", shortName="dontUseSoftClippedBases", required=false)
  var dontUseSoftClippedBases = false

  @Argument(doc="stand call conf.",
    fullName="stand_call_conf", shortName="stand_call_conf", required=false)
  var stand_call_conf: Double=20.0

  @Argument(doc="min mapping_quality score.",
    fullName="min_mapping_quality_score", shortName="min_mapping_quality_score", required=false)
  var min_mapping_quality_score: Int=20
  
  @Argument(doc="How many ways to scatter/gather",
    fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1
  
  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=false)
  var indels: Seq[File] = Seq() 
  

  def script() {
    var readyBams: Seq[File] = Seq()
    for (bam <- input) {
        if (qscript.nContigs < 0)
           qscript.nContigs = QScriptUtils.getNumberOfContigs(bam)
        //readyBams :+= run_printReads(run_bgsr(splitRNASeq(run_dedup(run_target_realign(bam)))))
        readyBams :+= run_printReads(run_bgsr(run_target_realign(splitRNASeq(run_dedup(bam)))))
	//readyBams :+= bam
    }
    genoTyper(readyBams)
  }
 
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 8
    this.isIntermediate = false
    
  }
  
  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
    //if(qscript.change255to60)
    //  this.read_filter :+= "ReassignOneMappingQuality"
    //if(allow_N_cigar_reads)
    //  this.unsafe = ValidationExclusion.TYPE.ALLOW_N_CIGAR_READS
  }
   
   def genoTyper(bams: Seq[File]) ={
    var uvcfFile = "UnifiedGenotyper.variant.vcf"
    var hvcfFile = "HaplotypeCaller.variant.vcf"
    
    val uvcfFileR = swapExt(uvcfFile, ".vcf", ".recal")
    val hvcfFileR = swapExt(hvcfFile, ".vcf", ".recal")
    val uvcfFileT = swapExt(uvcfFile, ".vcf", ".tr")
    val hvcfFileT = swapExt(hvcfFile, ".vcf", ".tr")
    val uvcfFileP = swapExt(uvcfFile, ".vcf", ".plot")
    val hvcfFileP = swapExt(hvcfFile, ".vcf", ".plot")
    if (qscript.caller == "UnifiedGenotyper") {
       add(unified_genotype(bams, uvcfFile))
      //add(vari_recal(uvcfFile,uvcfFileR,uvcfFileT,uvcfFileP))
    } else if (qscript.caller == "HaplotypeCaller") {
       add(haplo_genotype(bams, hvcfFile))
       // add(vari_recal(hvcfFile,hvcfFileR,hvcfFileT,hvcfFileP))
    } else if (qscript.caller == "UnifiedGenotyperAndHaplotypeCaller") {
       add(unified_genotype(bams, uvcfFile))
       add(haplo_genotype(bams, hvcfFile))
       // add(vari_recal(uvcfFile,uvcfFileR,uvcfFileT,uvcfFileP))
       // add(vari_recal(hvcfFile,hvcfFileR,hvcfFileT,hvcfFileP))
    }
  }  
  
  // java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar  -T VariantRecalibrator  -nt $NP  -R $R   -input my.variants.vcf  -recalFile my.snps.recal   -tranchesFile my.snps.tranches  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATKDATA/hapmap_3.3.b37.vcf  -resource:omni,known=false,training=true,truth=true,prior=12.0 $GATKDATA/1000G_omni2.5.b37.vcf  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATKDATA/1000G_phase1.snps.high_confidence.b37.vcf     -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp      -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0  -mode SNP  --maxGaussians 4  -rscriptFile recal.plots.R
  
//  java -Xmx4g -jar GenomeAnalysisTK.jar \
// *   -T VariantRecalibrator \
// *   -R reference/human_g1k_v37.fasta \
// *   -input NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf \
// *   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
// *   -resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
// *   -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 dbsnp_135.b37.vcf \
// *   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an InbreedingCoeff \
// *   -mode SNP \
// *   -recalFile path/to/output.recal \
// *   -tranchesFile path/to/output.tranches \
// *   -rscriptFile path/to/output.plots.R
 
   // didn't work. Not sure how to pass resource information to the class
   case class vari_recal (inVCF: File, recaFile: File, tranFile: File, plotFile:File) extends VariantRecalibrator with CommandLineGATKArgs{
    this.input_file = Seq(inVCF); 
    this.recalFile = recaFile
    this.tranchesFile = tranFile
    this.resource=qscript.dbSNP
    this.an=Seq("QD","MQRankSum","ReadPosRankSum","FS","MQ")
    this.tranche=Seq(100.0, 99.9, 90.0)
    //this.mod=Mod.SNP    // Error: value mod is not a member of gatkpipeline.this.vari_recal
    this.maxGaussians=4
    //this.rsscriptFile=plotFile
    //this.scatterCount = qscript.nContigs
    this.analysisName = inVCF + ".vari_recal"
    this.jobName = inVCF + ".vari_recal" 
  }
 
// java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar  -R $R  -T HaplotypeCaller  $bams  -recoverDanglingHeads  -dontUseSoftClippedBases  -stand_call_conf 20.0 -stand_emit_conf 20.0  --dbsnp $dbsnp  --output_mode EMIT_ALL_SITES  --genotyping_mode DISCOVERY  --emitRefConfidence BP_RESOLUTION  --out my.variants.vcf $interval"
     
  case class haplo_genotype ( inBams: Seq[File], vcfFile: File) extends HaplotypeCaller with CommandLineGATKArgs{
    this.input_file = inBams; 
    this.out = vcfFile
    this.stand_call_conf = qscript.stand_call_conf
    //this.stand_emit_conf = 20.0
    this.dontUseSoftClippedBases=qscript.dontUseSoftClippedBases
    this.genotyping_mode=GenotypingOutputMode.DISCOVERY
    this.min_mapping_quality_score=qscript.min_mapping_quality_score
    this.scatterCount = qscript.nContigs
    this.analysisName = "haplo_genotype"
    this.jobName =  "haplo_genotype"
  }
    // java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -glm BOTH -R $R $bams   -o my.variants.vcf -D $dbsnp  -stand_call_conf 30.0  -stand_emit_conf 10.0 $interval
  case class unified_genotype ( inBams: Seq[File], vcfFile: File) extends UnifiedGenotyper with CommandLineGATKArgs{
    this.input_file = inBams; 
    this.out = vcfFile //"my_unifiedtype.vcf"
    this.genotype_likelihoods_model=GenotypeLikelihoodsCalculationModel.Model.BOTH
    //this.stand_emit_conf = 10.0
    this.stand_call_conf = 30.0
    this.scatterCount = qscript.nContigs
    this.analysisName = "unified_genotype"
    this.jobName =  "unified_genotype"
  }
  
  def run_printReads (inBam: File) : File = {
    val recalTable = swapExt(inBam, ".bam", ".recal.table")
    val recalBam = swapExt(inBam, ".bam", ".recal.bam")
    add(printReads(inBam, recalTable,recalBam))
    recalBam
  }
   
   case class printReads (inBam: File, inRecalFile: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.BQSR = inRecalFile
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.out = outBam
    this.scatterCount = qscript.nContigs
    this.isIntermediate = false
    this.analysisName = outBam 
    this.jobName = outBam
  }
  
  // base quality recalibration   
  def run_bgsr(bam: File) : File = {
      val recalTable = swapExt(bam, ".bam", ".recal.table")
      add(bqsr(bam, recalTable))
      bam
  }
     
  case class bqsr ( inBam: File,  outRecalFile: File) extends BaseRecalibrator  with CommandLineGATKArgs {
    this.knownSites ++= qscript.dbSNP
    this.covariate ++= Seq("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "ContextCovariate")
    this.input_file :+= inBam
    this.out = outRecalFile
    this.scatterCount = qscript.nContigs
    this.memoryLimit = 4
    this.residentLimit = 6
    this.residentRequest = 6
    this.analysisName =  outRecalFile
  }
  
  
   def splitRNASeq(bam: File) : File ={
    if (splitNCigarReads) {
      val splitReadsBam = swapExt(bam,".bam",".splitReads.bam")
      add(splitRNAseqReads(bam, splitReadsBam))
      splitReadsBam
    }
    else {
      bam
    }
  }
  
   case class splitRNAseqReads(inBam: File, outBam: File) extends SplitNCigarReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.out = outBam
    this.scatterCount = qscript.nContigs
    //this.fix_misencoded_quality_scores = true
    this.read_filter :+= "ReassignOneMappingQuality"
    //this.unsafe = ValidationExclusion.TYPE.ALLOW_N_CIGAR_READS
  }
  
    
   // mark duplicate reads in bam file   
   def run_dedup(bam: File) : File = {
      val dedupedBam = swapExt(bam, ".bam", ".dedup.bam")
      val metricsFile = swapExt(bam, ".bam", ".dedup.metrics")
      add(dedup(bam, dedupedBam, metricsFile))
      dedupedBam
  }
  
  case class dedup ( inBam: File,  outBam: File, metricsFile: File) extends MarkDuplicates with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outBam
    this.metrics = metricsFile
    this.analysisName = outBam
  }
    
  //find indel region and re-align them
  def run_target_realign(inBam: File) : File =  {
      // both target and realign need File list to work, The original idea is to have multiple files for each sample
      val inBams: Seq[File] = Seq(inBam)
      val outInterval = swapExt(inBam, ".bam", ".realign.intervals")
      val outBam = swapExt(inBam, ".bam", ".realign.bam")
      add(target(inBams, outInterval))
      add(realign(inBams, outInterval, outBam))
      outBam
  }
  
   case class realign ( inBams: Seq[File], tIntervals: File,  outBam: File) extends IndelRealigner with CommandLineGATKArgs{
    this.input_file = inBams
    this.targetIntervals = tIntervals
    this.out = outBam
    this.known ++= qscript.dbSNP
    if (qscript.indels != null)
      this.known ++= qscript.indels
    this.consensusDeterminationModel = ConsensusDeterminationModel.USE_READS
    this.compress = 0
    this.scatterCount = qscript.nContigs
    this.analysisName = outBam 
    this.filter_reads_with_N_cigar = qscript.filter_reads_with_N_cigar
	this.jobName =  outBam 
  }
  
  // need dbSNP and nContig
  case class target ( inBams: Seq[File],  outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs{
    this.input_file = inBams
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.known ++= qscript.dbSNP
    this.filter_reads_with_N_cigar = qscript.filter_reads_with_N_cigar
    this.scatterCount = qscript.nContigs
    this.analysisName = outIntervals
  }
}
