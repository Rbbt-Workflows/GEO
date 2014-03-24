require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "GEO"

class TestClass < Test::Unit::TestCase
  def test_differential
    diff = GEO.job(:differential, "TEST",
            :dataset => "GDS4513",
            :main => "disease state=relapse").run

    assert diff.length > 0
  end

  def test_up_regulated
    genes = GEO.job(:up_genes, "TEST",
            :dataset => "GDS4513",
            :main => "disease state=relapse").run

    assert genes.length > 0
  end


  def test_down_regulated
    genes = GEO.job(:down_genes, "TEST",
            :dataset => "GDS4513",
            :main => "disease state=relapse").run

    assert genes.length > 0
  end
end


