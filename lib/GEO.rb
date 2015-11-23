# encoding: utf-8
require 'rbbt-util'
require 'rbbt/GE'
require 'rbbt/sources/organism'
require 'rbbt/resource'
require 'yaml'

module GEO
  extend Resource
  self.pkgdir = "geo"
  self.subdir = "arrays"

  GEO.claim GEO.root, :rake, Rbbt.share.install.GEO.Rakefile.find(:lib)

  def self.esearch(query)
    url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
      :db => "gds",
      :term => query + " AND gds[filter]",
      :retmax => "5000",
      :usehistory => "y",
    }

    url = url << '?' << Misc.hash2GET_params(params)
    Log.debug "Query GEO: #{url}"
    Nokogiri::XML(Open.open(url))
  end

  def self.efetch(query)
    url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
      :db => "gds",
      :term => query + " AND gds[filter]",
      :retmax => "5000",
      :usehistory => "y",
    }

    url = url << '?' << Misc.hash2GET_params(params)
    Log.debug "Fetch GEO: #{url}"
    Nokogiri::XML(Open.open(url))
  end

  def self.esummary(query_key,web_env)
    url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/summary.fcgi"
    params = {
      :db => "gds",
      :query_key => query_key,
      :web_env => web_env,
    }

    url = url << '?' << Misc.hash2GET_params(params)
    Log.debug "Summary GEO: #{url}"
    Nokogiri::XML(Open.open(url))
  end

  def self.query(query)
    require 'nokogiri'

    doc = esearch(query)
    web_env = doc.css('WebEnv')
    datasets = doc.css('eSearchResult IdList Id').collect{|e| "GDS" << e.content }
    datasets
  end

  def self.comparison_name(field, condition, control)
    condition = condition * " AND " if Array === condition
    control = control * " AND " if Array === control
    [[field, condition] * ": ", [field, control] * ": "] * " => "
  end

  def self.parse_comparison_name(name)
    field1, condition1, field2, condition2 = name.match(/(.*): (.*?) => (.*?): (.*)/).values_at(1, 2, 3, 4)
    condition1 = condition1.split(/ AND /) if condition1 =~ / AND /
    condition2 = condition2.split(/ AND /) if condition2 =~ / AND /

    [field1, condition1, field2, condition2]
  end

  def self.platform_info(platform)
    YAML.load(self[platform]['info.yaml'].produce.read)
  end

  def self.dataset_info(dataset)
    YAML.load(self[dataset]['info.yaml'].produce.read)
  end

  def self.is_control?(value)
    value.to_s.downcase =~ /\bcontrol\b/ or
    value.to_s.downcase =~ /\bwild/ or
    value.to_s.downcase =~ /\bnone\b/ 
  end

  def self.control_samples(dataset)
    info = dataset_info(dataset)
    subsets = info[:subsets]

    control_samples = []
    subsets.each do |type, values|
      control_samples.concat values.select{|value,samples| is_control? value}.collect{|value,samples| samples.split(",")}.flatten 
    end

    control_samples
  end

  module SOFT

    GDS_URL="ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/GDS/#DATASET#.soft.gz"
    GPL_URL="ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_platform/#PLATFORM#/#PLATFORM#_family.soft.gz"
    GSE_URL="ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/#SERIES#/#SERIES#_family.soft.gz"

    GSE_INFO = {
      :DELIMITER        => "\\^PLATFORM",
      :title         => "!Series_title",
      :channel_count => "!Sample_channel_count",
      :value_type    => "!Series_value_type",
      :platform      => "!Series_platform_id",
      :description   => "!Series_summary*",      # Join with \n 
    }
 
    GSE_SAMPLE_INFO = {
      :DELIMITER        => "\\^SAMPLE",
      :title         => "!Sample_title",
      :accession         => "!Sample_geo_accession",
      :channel_count => "!Sample_channel_count",
    }
   
    GDS_INFO = {
      :DELIMITER        => "\\^SUBSET",
      :value_type       => "!dataset_value_type",
      :channel_count    => "!dataset_channel_count",
      :platform         => "!dataset_platform",
      :reference_series => "!dataset_reference_series",
      :description      => "!dataset_description"
    }

    GDS_SUBSET_INFO = {
      :DELIMITER        => "!subset_.*|!dataset_value_type",
      :description => "!subset_description",
      :samples     => "!subset_sample_id*",
      :type        => "!subset_type",
    }

    GPL_INFO = { 
      :DELIMITER     => "!platform_table_begin",
      :organism      => "!Platform_organism",
      :count         => "!Platform_data_row_count"
    }

    # When multiple matches select most common, unless join is choosen
    def self.find_field(header, field, join = false)
      md = header.match(/#{ Regexp.quote field }\s*=\s*(.*)/i)
      return nil if md.nil? or md.captures.empty?

      case join
      when false, nil
        counts = Hash.new(0)
        md.captures.sort_by{|v| counts[v] += 1}.first
      when true
        md.captures * "\n"
      else
        md.captures * join
      end
    end

    def self.get_info(header, info)
      result = {}

      info.each do |key, field|
        next if key == :DELIMITER
        if field =~ /(.*)\*(.*)(\*)?$/
          value = find_field(header, $1, $2.empty? ? true : $2)
          value = value.to_i.to_s == value ? value.to_i : value
          if $3
            result[key] = value.split(',')
          else
            result[key] = value
          end
        else
          value = find_field(header, field, false)
          value = value.to_i.to_s == value ? value.to_i : value
          result[key] = value
        end
      end

      if result.empty?
        nil
      else
        result
      end
    end

    def self.parse_header(stream, info)
      header = ""
      while line = stream.readline
        line = Misc.fixutf8 line
        header << line
        break if line =~ /^#{info[:DELIMITER]}/i
        raise "Delimiter not found" if stream.eof?
      end

      get_info(header, info)
    end

    def self.guess_id(organism, codes)
      num_codes = codes.length
      best = nil
      best_count = 0
      new_fields = []
      field_counts = {}
      TmpFile.with_file(codes.to_s) do |codefile|

        codes.all_fields.each_with_index do |field,i|
          values = CMD.cmd("cat #{ codefile }|cut -f #{ i + 1 }| tr '|' '\\n'|grep [[:alpha:]]|sort -u").read.split("\n").reject{|code| code.empty?}

          new_field, count =  Organism.guess_id(organism, values)
          new_field ||= field
          count ||= 0
          new_field = "UNKNOWN(#{new_field})" unless count > (num_codes > 20000 ? 20000 : num_codes).to_f * 0.2 and count > values.uniq.length * 0.5

          Log.debug "Original field: #{ field }. New: #{new_field}. Count: #{ count }/#{num_codes}/#{values.uniq.length}"
          new_fields << new_field

          field_counts[new_field] = count

          if count > best_count
            best = new_field
            best_count = count
          end

        end

      end

      field_counts.delete(new_fields.first)
      [best, new_fields, field_counts.sort_by{|field, counts| counts}.collect{|field, counts| field}.compact]
    end

    #{{{ GPL

    def self.GPL(platform, directory)
      FileUtils.mkdir_p directory unless File.exists? directory

      code_file = File.join(directory, 'codes') 
      info_file = File.join(directory, 'info.yaml') 

      # Fix platforms with the '.\d' extension (eg. NM_020527.1)
      stream = CMD.cmd('sed \'s/\.[[:digit:]]\+\(\t\|$\)/\1/g;s/ *\/\/[^\t]*//g\'', :in =>  Open.open(GPL_URL.gsub('#PLATFORM#', platform), :nocache => true), :pipe => true)

      info = parse_header(stream, GPL_INFO)

      info[:code_file]      = code_file
      info[:data_directory] = directory

      Log.medium "Producing code file for #{ platform }"
      codes = TSV.open stream, :fix => proc{|l| l =~ /^!platform_table_end/i ? nil : l}, :header_hash => "", :sep2 => /\s*[|,]\s*/
      Log.low "Original fields: #{codes.key_field} - #{codes.fields * ", "}"

      best_field, all_new_fields, order = guess_id(Organism.default_code(Organism.organism(info[:organism])), codes)

      new_key_field, *new_fields = all_new_fields

      new_key_field = codes.key_field if new_key_field =~ /^UNKNOWN/

      codes.key_field = new_key_field.dup 
      codes.fields = new_fields.collect{|f| f.dup}

      Log.low "New fields: #{codes.key_field} - #{codes.fields * ", "}"

      Open.write(code_file, codes.reorder(:key, order).to_s(:sort))
      Open.write(info_file, info.to_yaml)

      info
    end

    def self.dataset_subsets(stream)
      text = ""
      while not (line = Misc.fixutf8(stream.gets)) =~ /!dataset_table_begin/
        text << line
      end

      subsets = text.split(/\^SUBSET/).collect do |chunk|
        get_info(chunk, GDS_SUBSET_INFO)
      end

      info = {}
      subsets.each do |subset|
        type = subset[:type]
        description = subset[:description]
        samples = subset[:samples]
        info[type] ||= {}
        info[type][description] = samples
      end

      info
    end

    def self.GDS(dataset, directory)
      FileUtils.mkdir_p directory unless File.exists? directory
      
      value_file = File.join(directory, 'values') 
      info_file = File.join(directory, 'info.yaml') 

      stream = Open.open(GDS_URL.gsub('#DATASET#', dataset), :nocache => true)

      info = parse_header(stream, GDS_INFO)
      info[:value_file]      = value_file
      info[:data_directory] = directory

      info[:subsets] = dataset_subsets(stream)
      platform = info[:platform]

      Log.medium "Producing values file for #{ dataset }"
      values = TSV.open stream, :fix => proc{|l| l =~ /^!dataset_table_end/i ? nil : l.gsub(/null/,'NA')}, :header_hash => "", :type => :list
      key_field = TSV.parse_header(GEO[platform].codes.produce.find).key_field
      values.key_field = key_field

      samples = values.fields.select{|f| f =~ /GSM/}

      Open.write(value_file, values.slice(samples).to_s(:sort))
      Open.write(info_file, info.to_yaml)

      info
    end

    def self.series_samples(stream)
      text = Misc.fixutf8(stream.read)

      values = nil

      sample_info = {}
      
      samples = []
      text.split(/\^SAMPLE/).each do |chunk|
        info = get_info(chunk, GSE_SAMPLE_INFO)
        sample = info[:accession]
        next if sample.nil?

        samples << sample

        chunk = chunk.encode "UTF-8"
        chunk = Misc.fixutf8 chunk
        sample_values = TSV.open(StringIO.new(chunk.match(/!sample_table_begin\n(.*)\n!sample_table_end/mi)[1].strip), :type => :list, :header_hash => '',:unnamed => true, :fields => [1])

        sample_values.fields = [sample]

        if values.nil?
          values = sample_values
        else
          values.attach sample_values
        end
        sample_info[sample] = info
      end

      [values, sample_info]
    end

    def self.GSE(series, directory)
      FileUtils.mkdir_p directory unless File.exists? directory

      value_file = File.join(directory, 'values') 
      info_file = File.join(directory, 'info.yaml') 

      stream = Open.open(GSE_URL.gsub('#SERIES#', series), :nocache => true)

      info = parse_header(stream, GSE_INFO)
      info[:value_file]      = value_file
      info[:data_directory] = directory

      Log.medium "Producing values file for #{ series }"
      values, sample_info = series_samples(stream)

      key_field = TSV.parse_header(GEO[info[:platform]]['codes'].open).key_field
      values.key_field = key_field

      info[:sample_info] ||= sample_info
      info[:channel_count] ||= sample_info.values.first[:channel_count]
      info[:value_type] ||= sample_info.values.first[:value_type]


      txt = values.dumper_stream.read
      Open.write(value_file, txt)
      Open.write(info_file, info.to_yaml)

      info
    end
  end



  def self.compare(dataset, field, condition, control, path)
    dataset_info = GEO[dataset]["info.yaml"].yaml

    platform = dataset_info[:platform]
    platform_info = GEO[platform]["info.yaml"].yaml

    log2       = ["count"].include? dataset_info[:value_type]
    samples    = dataset_info[:subsets]
    value_file = GEO[dataset].values.find.produce
    format     = GEO[platform].codes.open{|io| TSV.parse_header(io).key_field } 

    if Array === condition
      condition_samples = condition.collect{|cond| samples[field][cond].split ","}.flatten
    else
      condition_samples = samples[field][condition].split ","
    end

    if Array === control
      control_samples = control.collect{|cond| samples[field][cond].split ","}.flatten
    else
      control_samples = samples[field][control].split ","
    end

    GE.analyze(value_file, condition_samples, control_samples, log2, path, format)
  end

  def self.comparisons(all_conditions, controls)
    rest = all_conditions - controls
    comparisons = []

    if controls.any?
      controls.each do |control|
        rest.each do |ccase|
          comparisons << [ccase, control]
        end
      end
    else
      rest.each do |control|
        rest.each do |ccase|
          next if ccase == control
          comparisons << [ccase, control]
        end
      end
    end
    comparisons
  end

  def self.subset_comparisons(subsets)
    subset_comparisons = {}
    subsets.each do |subset,values|
      all_conditions = values.keys
      controls = all_conditions.select do |field| GEO.is_control? field end
      comparisons = GEO.comparisons all_conditions, controls

      subset_comparisons[subset] = {}
      comparisons.each do |ccase,control|
        case_samples = values[ccase]
        control_samples = values[control]
        case_samples = case_samples.split(",") if String === case_samples
        control_samples = control_samples.split(",") if String === control_samples
        next unless case_samples.length > 2 or control_samples.length > 2
        subset_comparisons[subset][[ccase,control]] = [case_samples, control_samples]
      end
    end

    subset_comparisons
  end

  def self.platform_key_field(platform)
    TSV.parse_header(GEO[platform].codes.produce.find).key_field
  end

  def self.dataset_key_field(database)
    platform_key_field dataset_info(database)[:platform]
  end

  def self.dataset_organism_code(dataset)
    info = dataset_info dataset
    Organism.default_code(Organism.organism(info[:organism]))
  end
  
  def self.dataset_comparisons(dataset)
    subset_comparisons dataset_info(dataset)[:subsets]
  end
end

