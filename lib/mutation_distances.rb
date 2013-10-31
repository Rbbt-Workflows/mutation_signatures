require 'rbbt/workflow'
require File.join(File.dirname(__FILE__), 'exome_positions')

module MutationSignatures
  input :mutations, :array, "Genomic Mutations"
  input :organism, :string, "Organism code"
  input :exome_only, :boolean, "Consider only exome bases", true
  task :mutation_distance => :tsv do |mutations,organism,exome_only|
    mutations = MutationSignatures.exome_position_at_genomic_positions(organism, mutations) if exome_only
    chr_mutations = {}
    mutations.each do |mutation|
      chr = mutation.split(":").first
      chr_mutations[chr] ||= []
      chr_mutations[chr] << mutation
    end
    chr_distances = TSV.setup({}, :key_field => "Chromosome Name", :fields => ["Distance"], :type => :flat)

    chr_mutations.each do |chr, mutations|
      last = 0
      chr_distances[chr] = []
      mutations.sort_by{|m| m.split(":")[1].to_i}.collect do |mutation|
        position = mutation.split(":")[1].to_i
        distance = position - last
        last = position
        chr_distances[chr] << [mutation, distance] * "@"
      end
    end

    chr_distances
  end

  dep :mutation_distance
  dep do |jobname, inputs| MutationSignatures.job(:exome_sizes, inputs[:organism], :organism => inputs[:organism]) end
  task :rainfall_plot => :binary do
    mutation_distance = step(:mutation_distance).load
    exome_only = step(:mutation_distance).info[:inputs][:exome_only]
    organism = step(:mutation_distance).info[:inputs][:organism]

    if exome_only
      exome_sizes = step(:exome_sizes).load
    else
      exome_sizes = Organism.chromosomes(organism).tsv(:type => :list, :fields => []).add_field("Chromosome Size"){|chromosome| begin File.size(Organism[organism]["chromosome_" << chromosome].find) rescue nil end }
    end

    FileUtils.mkdir files_dir unless File.exists? files_dir

    TmpFile.with_file(exome_sizes.to_s) do |sizes|
      TmpFile.with_file(mutation_distance.to_s) do |filename|
        script =<<-EOF
source('#{Rbbt.share.R["plots.R"].find(:lib)}')
data = rbbt.flat.tsv('#{filename}', as.is=T)
chr_sizes = rbbt.tsv('#{sizes}')
p <- rainfall_plot(data, chr_sizes)
ggsave('#{file('rainfall.png')}', p)
        EOF
        R.run script
      end
    end

    Open.read(file('rainfall.png'), :mode => "rb")
  end
  export_asynchronous :rainfall_plot
end
