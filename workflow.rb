require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/GE'
require 'GEO'

require 'rbbt/matrix'
require 'rbbt/matrix/differential'
require 'rbbt/matrix/barcode'
require 'rbbt/statistics/random_walk'

module GEO
  extend Workflow

  helper :organism do
    Organism.default_code("Hsa")
  end

  helper :dataset_info do |dataset|
    GEO[dataset]["info.yaml"].yaml
  end

  helper :identifiers do |dataset|
    if dataset =~ /^GDS|GSE/
      platform = dataset_info(dataset)[:platform]
      GEO[platform].codes
    else
      Organism.identifiers(organism)
    end
  end

  helper :translate do |dataset,keys,format|
    log :translate, "Translating keys to #{ format }"

    if dataset =~ /^GDS|GSE/
      platform = dataset_info(dataset)[:platform]
      codes = GEO[platform].codes
      if TSV.parse_header(codes).all_fields.include? format
        index = codes.index :target => format
      end
    end

    index ||= Organism.identifiers(organism).index :target => format, :persist => true

    keys.collect{|k| index[k] }
  end

  helper :matrix do |dataset|
    if dataset =~ /^GDS|GSE/
      dataset_info = dataset_info dataset

      samples       = dataset_info[:subsets]
      value_type    = dataset_info[:value_type]

      platform = dataset_info[:platform]
      format     = GEO[platform].codes.open{|io| TSV.parse_header(io).key_field } 

      log :geo, "Processing #{ dataset } files from GEO" unless File.exists?(GEO[dataset].values.find)
      data = GEO[dataset].values 

      log :matrix, "Setting matrix for #{ dataset }"
      matrix = Matrix.new data.find.produce, samples, value_type, format
      matrix.subsets = dataset_info(dataset)[:subsets]
      matrix
    else
      raise "Matrix not found: #{ dataset }" unless File.exists?(dataset)
      format = TSV.parse_header(dataset).key_field

      log :matrix, "Setting matrix for #{ dataset }"
      matrix = Matrix.new Path.setup(dataset.dup), nil, nil, format
      matrix.subsets = dataset_info(dataset)[:subsets]
      matrix
    end
  end

  helper :matrix_to_gene do |matrix, dataset|
    log :translating, "Translating #{ dataset } matrix to known gene ids"
    begin
      subsets = matrix.subsets
      matrix.to_gene(identifiers(dataset)) 
    rescue FieldNotFoundError
      raise ParameterException, "Could not identify probes in #{ dataset }. Cannot translate to gene -- #{identifiers(dataset).find}"
    end 
  end


  input :dataset, :string, "Dataset code", nil
  task :dataset_info => :yaml do |dataset|
    raise ParameterException, :dataset if dataset.nil?
    GEO.dataset_info(dataset)
  end

  input :dataset, :string, "Dataset code", nil
  task :platform => :string do |dataset|
    GEO.dataset_info(dataset)[:platform]
  end

  input :platform, :string, "Platform code", nil
  task :platform_info => :yaml do |platform|
    GEO.platform_info(platform)
  end



  input :query, :string, "Query"
  task :query => :array do |query|
    raise ParameterException, "No query" if query.nil? or query.empty?
    GEO.query query
  end
  export_asynchronous :query
  
  input :dataset, :string, "Dataset code", nil
  task :sample_info => :tsv do |dataset|
    dataset_info = dataset_info dataset

    if dataset_info.include? :sample_info
      sample_info = dataset_info[:sample_info]

      tsv = TSV.setup({}, :key_field => "Sample", :fields => ["Title"], :type => :single)

      sample_info.each do |sample,info|
        tsv[sample] = info[:title]
      end

      tsv
    else
      subsets = dataset_info[:subsets]

      tsv = TSV.setup({}, :key_field => "Sample", :fields => [], :type => :list)
      subsets.each do |factor,values|
        factor_samples = {}
        values.each do |value,samples|
          samples.split(",").each do |sample|
            factor_samples[sample] = value
          end
        end
        factor_samples.each do |sample,value|
          tsv[sample] ||= []
          tsv[sample] << value
        end
        tsv.fields << factor
      end

      tsv
    end
  end
  export_asynchronous :sample_info

  input :dataset, :string, "Dataset code", nil
  input :to_gene, :boolean, "Average probes by gene", false
  task :matrix => :tsv do |dataset,to_gene|
    raise ParameterException, "No dataset provided" if dataset.nil?
    matrix = matrix dataset

    matrix = matrix_to_gene matrix, dataset if to_gene

    FileUtils.cp matrix.data_file, path
  end
  export_asynchronous :matrix

  input :dataset, :string, "Dataset code", nil
  input :main, :string, "Main samples", nil
  input :contrast, :string, "Contrat samples", nil
  input :to_gene, :boolean, "Average probes by gene", false
  task :differential => :tsv do |dataset,main, contrast,to_gene|
    raise ParameterException, "No dataset provided" if dataset.nil?
    matrix = matrix dataset

    if main and main[0] == "/"
      info = dataset_info(dataset) 
      re = Regexp.new(main[1..-2])
      main = info[:sample_info].select{|k,v| v[:title] =~ re}.collect{|k,v| k}
    end

    if contrast and contrast[0] == "/"
      info = dataset_info(dataset) 
      re = Regexp.new(contrast[1..-2])
      contrast = info[:sample_info].select{|k,v| v[:title] =~ re}.collect{|k,v| k}
    end

    matrix = matrix_to_gene matrix, dataset if to_gene

    log :differential, "Running differential for #{ dataset } -- #{Misc.fingerprint [main, contrast]}"
    matrix.differential main, contrast, file(:result)
    FileUtils.mv file(:result), path
    nil
  end
  export_asynchronous :differential

  dep :differential
  input :threshold, :float, "Cuttoff", 0.005
  returns "Ensembl Gene ID"
  task :up_genes => :array do |threshold|
    diff = step(:differential).load
    dataset = step(:differential).info[:inputs][:dataset]
    to_gene = step(:differential).info[:inputs][:to_gene]

    diff.unnamed = true
    keys = diff.select("adjusted.p.values"){|p| p > 0 and p < threshold }.keys

    to_gene ? translate(dataset, keys, "Ensembl Gene ID") : keys
  end
  export_asynchronous :up_genes

  dep :differential
  input :threshold, :float, "Cuttoff", 0.005
  returns "Ensembl Gene ID"
  task :down_genes => :array do |threshold|
    diff = step(:differential).load
    dataset = step(:differential).info[:inputs][:dataset]
    to_gene = step(:differential).info[:inputs][:to_gene]

    diff.unnamed = true
    keys = diff.select("adjusted.p.values"){|p| p < 0 and p.abs < threshold }.keys

    to_gene ? translate(dataset, keys, "Ensembl Gene ID") : keys
  end
  export_asynchronous :down_genes


  input :dataset, :string, "Dataset code", nil
  input :to_gene, :boolean, "Average probes by gene", false
  input :standard_deviations, :float, "Standard deviations from the lower mode", 2
  task :barcode => :tsv do |dataset,to_gene, standard_deviations|
    raise ParameterException, "No dataset provided" if dataset.nil?
    matrix = matrix dataset

    matrix = matrix_to_gene matrix, dataset if to_gene

    matrix.barcode path, standard_deviations
  end
  export_asynchronous :barcode

  input :dataset, :string, "Dataset"
  input :to_gene, :boolean, "Transform to gene", true
  task :signatures => :tsv do |dataset,to_gene|
    subset_comparisons = GEO.dataset_comparisons dataset
    comparison_jobs = {}
    subset_comparisons.each do |subset,values|
      values.each do |comparison,samples|
        main, contrast = comparison
        main_samples, contrast_sammples = samples
        comparison_name = subset + ": " + [main, contrast] * " => "
        comparison_jobs[comparison_name] = GEO.job(:differential, comparison_name, :dataset => dataset, :main => main_samples, :contrast => contrast_sammples, :to_gene => to_gene)
      end
    end

    Misc.bootstrap comparison_jobs.values do |job|
      job.produce
    end


    if to_gene
      probe_id = "Ensembl Gene ID"
    else
      probe_id = GEO.dataset_key_field dataset
    end

    index_file = file('index')
    database = Persist.open_tokyocabinet index_file, true, :flat, "HDB"
    TSV.setup(database, :key_field => "Signature", :fields => [probe_id], :type => :flat, :unnamed => true, :namepace => GEO.dataset_organism_code(dataset))
    comparison_jobs.each do |name,job|
      database[name] = job.load.sort_by("t.values").collect{|k,v| k}
    end

    database.read

    database.to_s
  end

  dep GEO, :signatures do |jobname,options|
    GEO.job(:signatures, "Basic Signatures", options)
  end
  input :dataset, :string, "Dataset"
  input :up_genes, :array, "Up genes"
  input :down_genes, :array, "Down genes"
  input :permutations, :integer, "Number of permutation", 100_000
  task :rank_query => :tsv do |dataset,up_genes,down_genes,permutations|
    signature_job = step(:signatures)
    to_gene = signature_job.inputs[:to_gene]

    if to_gene
      index = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :persist  => true, :order => true, :unnamed => true
      up = index.chunked_values_at(up_genes) if up_genes
      down = index.chunked_values_at(down_genes) if down_genes
    else
      up = up_genes
      down = down_genes
    end

    index_file = signature_job.file('index')
    database = Persist.open_tokyocabinet index_file, false, :flat, "HDB"
    TSV.setup(database)

    dumper = TSV::Dumper.new :key_field => "Signature", :fields => ["P-value", "Hits"], :type => :double
    dumper.init
    FileUtils.mkdir_p files_dir
    TSV.traverse database, :cpus => 1, :into => dumper, :bar => true do |signature, list|
      list = OrderedList.setup(list)
      img_width ||= list.length / 100
      pvalue = list.pvalue(up, 0.2, :persist_permutations => true, :permutations => permutations)
      matches = (list & up)
      list.draw_hits(up, file(signature + '.png'), :width => img_width)
      [signature, [pvalue, matches]]
    end
  end



  dep GEO, :rank_query do |jobname, options|
    jobs = options[:datasets].collect do |dataset|
      GEO.job(:rank_query, jobname + ": " + dataset, options.merge(:dataset => dataset))
    end
    Misc.bootstrap jobs, 10 do |job|
      job.produce
    end
    jobs
  end
  input :datasets, :array, "List of GDS"
  task :rank_query_batch => :tsv do
    res = TSV.setup({}, :key_field => "Signature", :fields => ["P-value", "Hits"], :type => :single, :unnamed => true)

    dependencies.each do |job|
      partial_res = job.load
      dataset = job.inputs[:dataset]
      partial_res.each do |signature, value|
        signature = dataset + ": " + signature
        res[signature] = value
      end
    end

    res
  end
end
