require 'rbbt/workflow'
require 'rbbt/fix_width_table'
require 'rbbt/sources/organism'

module MutationSignatures
  def self.exon_ranges_for_chr(organism)
    @@exon_range_tsv ||= {}
    @@exon_range_tsv[organism] ||= Organism.exons(organism).tsv :persist => true, :fields => ["Exon Chr Start", "Exon Chr End", "Chromosome Name"], :type => :list, :unnamed => true
  end

  def self.exon_offsets_in_exome_for_chr(organism, chr)
    @@exon_offsets_in_exome_for_chr ||= {}
    @@exon_offsets_in_exome_for_chr[organism] ||= {}
    if @@exon_offsets_in_exome_for_chr[organism][chr].nil?
      fwt = Persist.persist("Exome offset of exons in #{ chr }[#{ organism }]", :fwt) do
        exon_ranges = exon_ranges_for_chr(organism)
        gaps = 0
        last = nil
        value_size = 0

        exome_positions = []
        exon_ranges.filter
        exon_ranges.add_filter("field:Chromosome Name", chr)
        exon_ranges.sort_by("Exon Chr Start"){|k,v| v.to_i}.each do |exon, range|
          start, eend = range.collect{|v| v.to_i}
          case
          when last.nil?
            gaps += start
          when last < start
            gaps += start - last - 1
          end

          last = eend if last.nil? or eend > last

          value = gaps.to_s
          value_size = value.length if value.length > value_size
          exome_positions << [value, start]
        end
        exon_ranges.pop_filter
        exon_ranges.reset_filters
        fwt = FixWidthTable.new :memory, value_size, false
        fwt.add_point exome_positions

        fwt
      end
      @@exon_offsets_in_exome_for_chr[organism][chr] = fwt
    end
    @@exon_offsets_in_exome_for_chr[organism][chr]
  end

  def self.exome_size(organism, chr)
    exon_ranges = exon_ranges_for_chr(organism)
    exome_positions = []
    exon_ranges.filter
    exon_ranges.add_filter("field:Chromosome Name", chr)
    begin
      eend = exon_ranges.sort_by("Exon Chr End"){|k,v| v.to_i}.last.last[1].to_i
    rescue
      Log.debug("Error calculating exome sizes for #{ organism } #{ chr }")
    end

    fwt = exon_offsets_in_exome_for_chr(organism, chr)
    eend_pos = fwt.pos(fwt.size - 1).to_s
    eend_gap = fwt.value(fwt.size - 1).to_s

    eend_pos.to_i - eend_gap.to_i
  end

  input :organism, :string, "Organism code", "Hsa"
  def self.exome_sizes(organism)
    tsv = TSV.setup({}, :key_field => "Chromosome Name", :fields => ["Exome Size"], :type => :single, :namespace => organism, :cast => :to_i)
    Organism.chromosomes(organism).tsv.keys.each do |chromosome|
      tsv[chromosome] = exome_size(organism, chromosome)
    end
    tsv
  end
  task :exome_sizes => :tsv

  desc "Exome position for chromosome position"
  input :organism, :string, "Organism code", "Hsa"
  input :chromosome, :string, "Chromosome name"
  input :positions, :array, "Positions"
  def self.exome_position_at_chr_positions(organism, chromosome, positions)
    fwt = exon_offsets_in_exome_for_chr(organism, chromosome)
    begin
      positions.collect{|p|
        p = p.to_i
        closest = fwt.closest(p)
        gap = fwt.value(closest).to_i
        p - gap
      }
    rescue
      raise $!
    end
  end
  task :exome_position_at_chr_positions => :array
  export_exec :exome_position_at_chr_positions

  desc "Exome position for genomic position"
  input :organism, :string, "Organism code", "Hsa"
  input :positions, :array, "Positions Chr:Position (e.g. 19:54646887). Separator can be ':', space or tab. Extra fields are ignored"
  def self.exome_position_at_genomic_positions(organism, positions)
    chr_positions = {}

    log :parsing, "Parsing positions"
    positions.each do |position|
      chr, pos = position.split(/[\s:\t]/)
      chr.sub!(/chr/,'')
      chr_positions[chr] ||= []
      chr_positions[chr] << pos
    end

    log :processing, "Processing chromosome positions"
    chr_exome_positions = {}
    chr_positions.each do |chr, list|
      chr_exome_positions[chr] = exome_position_at_chr_positions(organism, chr, list)
    end

    log :loading, "Loading results"
    positions.collect do |position|
      chr, pos, *rest = position.split(/[\s:\t]/)
      chr.sub!(/chr/,'')
      exome_position = chr_exome_positions[chr].shift
      ([chr, exome_position].concat rest) * ":"
    end
  end
  task :exome_position_at_genomic_positions=> :array
  export_synchronous :exome_position_at_genomic_positions
end
