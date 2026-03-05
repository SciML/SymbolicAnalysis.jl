<script setup lang="ts">
import { ref, computed, onMounted, onUnmounted } from 'vue'
import MarkdownIt from 'markdown-it'
import mathjax3 from 'markdown-it-mathjax3'

const md = MarkdownIt({ html: true }).use(mathjax3)

/* ------------------------------------------------------------------ */
/*  Built-in data & column configs (no external props needed)          */
/* ------------------------------------------------------------------ */

import dcpAtoms from '../../../data/dcp_atoms.json'
import dcpPowerAtoms from '../../../data/dcp_power_atoms.json'
import dgcpSpdAtoms from '../../../data/dgcp_spd_atoms.json'
import dgcpLorentzAtoms from '../../../data/dgcp_lorentz_atoms.json'

interface ColumnDef {
  key: string
  label: string
  filterable?: boolean
  sortable?: boolean
  icon?: boolean
  monoIcon?: boolean
  math?: boolean
}

const dataSets: Record<string, Record<string, string>[]> = {
  dcp: dcpAtoms as Record<string, string>[],
  'dcp-power': dcpPowerAtoms as Record<string, string>[],
  'dgcp-spd': dgcpSpdAtoms as Record<string, string>[],
  'dgcp-lorentz': dgcpLorentzAtoms as Record<string, string>[],
}

const columnConfigs: Record<string, ColumnDef[]> = {
  dcp: [
    { key: 'atom',         label: 'Atom',         sortable: true,  filterable: false },
    { key: 'meaning',      label: 'Meaning',      sortable: false, filterable: false, math: true },
    { key: 'domain',       label: 'Domain',       sortable: false, filterable: false, math: true },
    { key: 'curvature',    label: 'Curvature',    sortable: true,  filterable: true, icon: true },
    { key: 'monotonicity', label: 'Monotonicity', sortable: true,  filterable: true, monoIcon: true },
  ],
  'dcp-power': [
    { key: 'condition',    label: 'Condition',    sortable: false, filterable: false, math: true },
    { key: 'meaning',      label: 'Meaning',      sortable: false, filterable: false, math: true },
    { key: 'domain',       label: 'Domain',       sortable: false, filterable: false, math: true },
    { key: 'curvature',    label: 'Curvature',    sortable: true,  filterable: true, icon: true },
    { key: 'monotonicity', label: 'Monotonicity', sortable: true,  filterable: true, monoIcon: true },
  ],
  'dgcp-spd': [
    { key: 'atom',         label: 'Atom',         sortable: true,  filterable: false },
    { key: 'meaning',      label: 'Meaning',      sortable: false, filterable: false, math: true },
    { key: 'curvature',    label: 'G-Curvature',  sortable: true,  filterable: true, icon: true },
    { key: 'monotonicity', label: 'Monotonicity', sortable: true,  filterable: true, monoIcon: true },
  ],
  'dgcp-lorentz': [
    { key: 'atom',         label: 'Atom',         sortable: true,  filterable: false },
    { key: 'meaning',      label: 'Meaning',      sortable: false, filterable: false, math: true },
    { key: 'curvature',    label: 'G-Curvature',  sortable: true,  filterable: true, icon: true },
    { key: 'monotonicity', label: 'Monotonicity', sortable: true,  filterable: true, monoIcon: true },
  ],
}

/* ------------------------------------------------------------------ */
/*  Props — only need a source key string                              */
/* ------------------------------------------------------------------ */

const props = defineProps({
  src: { type: String, required: true },
})

const data = computed(() => dataSets[props.src] ?? [])
const columns = computed(() => columnConfigs[props.src] ?? [])

/* ------------------------------------------------------------------ */
/*  LaTeX rendering via MathJax (markdown-it-mathjax3)                 */
/* ------------------------------------------------------------------ */

function renderLatex(tex: string): string {
  if (!tex) return ''
  try {
    return md.render(`$${tex}$`).replace(/^<p>|<\/p>\s*$/g, '')
  } catch {
    return tex
  }
}

/* ------------------------------------------------------------------ */
/*  Docstring links — map atom name → functions page anchor            */
/* ------------------------------------------------------------------ */

const docLinks: Record<string, string> = {
  // DCP atoms
  'dotsort': 'SymbolicAnalysis.dotsort-Tuple{AbstractVector, AbstractVector}',
  'invprod': 'SymbolicAnalysis.invprod-Tuple{AbstractVector}',
  'eigsummax': 'SymbolicAnalysis.eigsummax-Tuple{LinearAlgebra.Symmetric, Int64}',
  'eigsummin': 'SymbolicAnalysis.eigsummin-Tuple{LinearAlgebra.Symmetric, Int64}',
  'matrix_frac': 'SymbolicAnalysis.matrix_frac-Tuple{AbstractVector, AbstractMatrix}',
  'perspective(f, x, s)': 'SymbolicAnalysis.perspective-Tuple{Function, Any, Real}',
  'quad_form': 'SymbolicAnalysis.quad_form-Tuple{AbstractVector, AbstractMatrix}',
  'quad_over_lin (array)': 'SymbolicAnalysis.quad_over_lin-Tuple{Real, Real}',
  'quad_over_lin (scalar)': 'SymbolicAnalysis.quad_over_lin-Tuple{Real, Real}',
  'sum_largest': 'SymbolicAnalysis.sum_largest-Tuple{AbstractMatrix, Integer}',
  'sum_smallest': 'SymbolicAnalysis.sum_smallest-Tuple{AbstractMatrix, Integer}',
  'trinv': 'SymbolicAnalysis.trinv-Tuple{AbstractMatrix}',
  'tv (vector)': 'SymbolicAnalysis.tv-Tuple{AbstractVector{<:Real}}',
  'tv (matrix)': 'SymbolicAnalysis.tv-Tuple{AbstractVector{<:AbstractMatrix}}',
  'huber': 'SymbolicAnalysis.huber',
  'lognormcdf': 'SymbolicAnalysis.lognormcdf-Tuple{Real}',
  // DGCP SPD atoms
  'conjugation': 'SymbolicAnalysis.conjugation-Tuple{Any, Any}',
  'scalar_mat': 'SymbolicAnalysis.scalar_mat',
  'sdivergence': 'SymbolicAnalysis.sdivergence-Tuple{Any, Any}',
  // quad_form already covered by DCP entry above
  'log_quad_form': 'SymbolicAnalysis.log_quad_form-Tuple{Vector{<:Number}, Matrix}',
  'schatten_norm': 'SymbolicAnalysis.schatten_norm',
  'sum_log_eigmax': 'SymbolicAnalysis.sum_log_eigmax-Tuple{Function, AbstractMatrix, Int64}',
  'affine_map': 'SymbolicAnalysis.affine_map-Tuple{typeofSymbolicAnalysis.conjugation, Matrix, Matrix, Matrix}',
  'hadamard_product': 'SymbolicAnalysis.hadamard_product-Tuple{AbstractMatrix, AbstractMatrix}',
  // DGCP Lorentz atoms
  'lorentz_log_barrier': 'SymbolicAnalysis.lorentz_log_barrier-Tuple{AbstractVector}',
  'lorentz_homogeneous_quadratic': 'SymbolicAnalysis.lorentz_homogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector}',
  'lorentz_homogeneous_diagonal': 'SymbolicAnalysis.lorentz_homogeneous_diagonal-Tuple{AbstractVector, AbstractVector}',
  'lorentz_least_squares': 'SymbolicAnalysis.lorentz_least_squares-Tuple{AbstractMatrix, AbstractVector, AbstractVector}',
  'lorentz_transform': 'SymbolicAnalysis.lorentz_transform-Tuple{AbstractMatrix, AbstractVector}',
}

function getDocLink(atomName: string): string | null {
  const anchor = docLinks[atomName]
  return anchor ? `./functions.html#${anchor}` : null
}

/* ------------------------------------------------------------------ */
/*  Search & filter state                                              */
/* ------------------------------------------------------------------ */

const searchQuery = ref('')

// Multi-select: each filterable column has a Set of selected values
// Empty set = show all (no filter active)
const selectedFilters = ref<Record<string, Set<string>>>({})

function toggleFilterValue(colKey: string, value: string) {
  if (!selectedFilters.value[colKey]) {
    selectedFilters.value[colKey] = new Set()
  }
  const s = selectedFilters.value[colKey]
  if (s.has(value)) {
    s.delete(value)
  } else {
    s.add(value)
  }
  // Force reactivity
  selectedFilters.value = { ...selectedFilters.value }
}

function clearFilter(colKey: string) {
  selectedFilters.value[colKey] = new Set()
  selectedFilters.value = { ...selectedFilters.value }
}

function isFilterActive(colKey: string): boolean {
  const s = selectedFilters.value[colKey]
  return !!s && s.size > 0
}

// Distinct options per filterable column (derived from data)
const filterOptions = computed(() => {
  const opts: Record<string, string[]> = {}
  for (const col of columns.value) {
    if (!col.filterable) continue
    const unique = [...new Set(data.value.map(r => r[col.key]))]
      .filter(Boolean)
      .sort()
    opts[col.key] = unique
  }
  return opts
})

/* ------------------------------------------------------------------ */
/*  Column header filter dropdown state                                */
/* ------------------------------------------------------------------ */

const openFilterKey = ref<string | null>(null)

function toggleFilterDropdown(colKey: string, event: MouseEvent) {
  event.stopPropagation()
  openFilterKey.value = openFilterKey.value === colKey ? null : colKey
}

function closeAllDropdowns() {
  openFilterKey.value = null
}

// Close dropdown when clicking outside
function handleDocumentClick(e: MouseEvent) {
  const target = e.target as HTMLElement
  if (!target.closest('.header-filter-dropdown') && !target.closest('.filter-toggle-btn')) {
    closeAllDropdowns()
  }
}

onMounted(() => document.addEventListener('click', handleDocumentClick))
onUnmounted(() => document.removeEventListener('click', handleDocumentClick))

/* ------------------------------------------------------------------ */
/*  Sort state                                                         */
/* ------------------------------------------------------------------ */

const sortKey = ref<string | null>(null)
const sortDir = ref<'asc' | 'desc'>('asc')

function toggleSort(key: string) {
  if (sortKey.value === key) {
    if (sortDir.value === 'asc') sortDir.value = 'desc'
    else { sortKey.value = null; sortDir.value = 'asc' }
  } else {
    sortKey.value = key
    sortDir.value = 'asc'
  }
}

function sortIndicator(key: string) {
  if (sortKey.value !== key) return '⇅'
  return sortDir.value === 'asc' ? '↑' : '↓'
}

/* ------------------------------------------------------------------ */
/*  Filtered + sorted rows                                             */
/* ------------------------------------------------------------------ */

const nameKey = computed(() => {
  const first = columns.value[0]
  return first ? first.key : 'atom'
})

const filteredRows = computed(() => {
  let rows = data.value

  // Text search
  if (searchQuery.value.trim()) {
    const q = searchQuery.value.trim().toLowerCase()
    rows = rows.filter(r => (r[nameKey.value] ?? '').toLowerCase().includes(q))
  }

  // Multi-select filters
  for (const [key, selected] of Object.entries(selectedFilters.value)) {
    if (selected.size > 0) {
      rows = rows.filter(r => selected.has(r[key]))
    }
  }

  // Sort
  if (sortKey.value) {
    const k = sortKey.value
    const dir = sortDir.value === 'asc' ? 1 : -1
    rows = [...rows].sort((a, b) => {
      const va = (a[k] ?? '').toLowerCase()
      const vb = (b[k] ?? '').toLowerCase()
      return va < vb ? -dir : va > vb ? dir : 0
    })
  }

  return rows
})

/* ------------------------------------------------------------------ */
/*  Curvature → icon mapping                                           */
/* ------------------------------------------------------------------ */

const curvatureIcons: Record<string, { color: string; darkColor: string; path: string; label: string }> = {
  convex:   { color: '#4063D8', darkColor: '#6B8BFF', path: 'M 4 2 Q 12 22 20 2',  label: 'Convex'   },
  concave:  { color: '#389826', darkColor: '#5BC848', path: 'M 4 20 Q 12 0 20 20', label: 'Concave'  },
  affine:   { color: '#CB3C33', darkColor: '#FF6B61', path: 'M 4 20 L 20 4',       label: 'Affine'   },
  gconvex:  { color: '#9558B2', darkColor: '#B87FD4', path: 'M 4 2 Q 12 22 20 2',  label: 'GConvex'  },
  glinear:  { color: '#9558B2', darkColor: '#B87FD4', path: 'M 4 20 L 20 4',       label: 'GLinear'  },
  gconcave: { color: '#9558B2', darkColor: '#B87FD4', path: 'M 4 20 Q 12 0 20 20', label: 'GConcave' },
}

function curvatureKey(raw: string): string | null {
  if (!raw) return null
  const l = raw.toLowerCase().replace(/\s+/g, '')
  if (l === 'gconvex') return 'gconvex'
  if (l === 'glinear') return 'glinear'
  if (l === 'gconcave') return 'gconcave'
  if (l === 'convex') return 'convex'
  if (l === 'concave') return 'concave'
  if (l === 'affine') return 'affine'
  return null
}

/* ------------------------------------------------------------------ */
/*  Monotonicity → icon mapping                                        */
/* ------------------------------------------------------------------ */

const monotonicityIcons: Record<string, { color: string; darkColor: string; path: string; arrowPath?: string; label: string; dashed: boolean }> = {
  increasing:             { color: '#389826', darkColor: '#5BC848', path: 'M 5 19 L 19 5',  arrowPath: 'M 13 5 L 19 5 L 19 11',    label: 'Increasing',     dashed: false },
  decreasing:             { color: '#CB3C33', darkColor: '#FF6B61', path: 'M 5 5 L 19 19',  arrowPath: 'M 13 19 L 19 19 L 19 13',  label: 'Decreasing',     dashed: false },
  anymono:                { color: '#888888', darkColor: '#AAAAAA', path: 'M 4 12 L 20 12',                                          label: 'Non-monotonic',  dashed: false },
  increasing_if_positive: { color: '#389826', darkColor: '#5BC848', path: 'M 5 19 L 19 5',  arrowPath: 'M 13 5 L 19 5 L 19 11',    label: 'Incr. if pos.',  dashed: true  },
  gincreasing:            { color: '#9558B2', darkColor: '#B87FD4', path: 'M 5 19 L 19 5',  arrowPath: 'M 13 5 L 19 5 L 19 11',    label: 'GIncreasing',    dashed: true  },
  gdecreasing:            { color: '#9558B2', darkColor: '#B87FD4', path: 'M 5 5 L 19 19',  arrowPath: 'M 13 19 L 19 19 L 19 13',  label: 'GDecreasing',    dashed: true  },
  ganymono:               { color: '#9558B2', darkColor: '#B87FD4', path: 'M 4 12 L 20 12',                                          label: 'GAnyMono',       dashed: true  },
}

function monoKey(raw: string): string | null {
  if (!raw) return null
  const l = raw.toLowerCase().replace(/\s+/g, '')
  if (l === 'increasing') return 'increasing'
  if (l === 'decreasing') return 'decreasing'
  if (l === 'anymono' || l === 'any') return 'anymono'
  if (l.startsWith('increasing_if_positive') || l.startsWith('increasingifpositive') || l.startsWith('incr.ifpos.')) return 'increasing_if_positive'
  if (l === 'gincreasing') return 'gincreasing'
  if (l === 'gdecreasing') return 'gdecreasing'
  if (l === 'ganymono' || l === 'gany') return 'ganymono'
  return null
}

interface MonoPart {
  key: string | null
  text: string
}

function parseMonotonicity(raw: string): MonoPart[] {
  if (!raw || raw === '—') return [{ key: null, text: '—' }]
  const tupleMatch = raw.match(/^\((.+)\)$/)
  if (tupleMatch) {
    const inner = tupleMatch[1]
    const parts: string[] = []
    let depth = 0, start = 0
    for (let i = 0; i < inner.length; i++) {
      if (inner[i] === '(') depth++
      else if (inner[i] === ')') depth--
      else if (inner[i] === ',' && depth === 0) {
        parts.push(inner.substring(start, i).trim())
        start = i + 1
      }
    }
    parts.push(inner.substring(start).trim())
    return parts.map(p => ({ key: monoKey(p), text: p }))
  }
  return [{ key: monoKey(raw), text: raw }]
}
</script>

<template>
  <div :class="['atom-table-wrapper', src]">
    <!-- Table -->
    <div class="atom-table-scroll">
      <table class="atom-table">
        <thead>
          <tr>
            <th
              v-for="col in columns"
              :key="col.key"
              :class="['col-' + col.key, { sortable: col.sortable !== false, 'has-filter': col.filterable }]"
            >
              <div class="th-content">
                <!-- Column label -->
                <span class="th-label">{{ col.label }}</span>

                <!-- Sort & filter controls row -->
                <div class="th-controls" v-if="col.sortable !== false || col.filterable">
                  <!-- Sort button -->
                  <button
                    v-if="col.sortable !== false"
                    class="sort-btn"
                    :title="`Sort by ${col.label}`"
                    @click="toggleSort(col.key)"
                  >
                    <svg viewBox="0 0 16 16" width="11" height="11" fill="currentColor">
                      <path d="M3 7l5-5 5 5H3z" :opacity="sortKey === col.key && sortDir === 'asc' ? 1 : 0.3"/>
                      <path d="M3 9l5 5 5-5H3z" :opacity="sortKey === col.key && sortDir === 'desc' ? 1 : 0.3"/>
                    </svg>
                  </button>

                  <!-- Filter toggle button -->
                  <button
                    v-if="col.filterable"
                    class="filter-toggle-btn"
                    :class="{ active: isFilterActive(col.key), open: openFilterKey === col.key }"
                    :title="`Filter by ${col.label}`"
                    @click="toggleFilterDropdown(col.key, $event)"
                  >
                    <svg viewBox="0 0 16 16" width="11" height="11" fill="currentColor">
                      <path d="M1 2h14l-5 6v5l-4 2V8L1 2z"/>
                    </svg>
                  </button>
                </div>
              </div>

              <!-- Dropdown panel (positioned absolutely below the header) -->
              <div
                v-if="col.filterable && openFilterKey === col.key"
                class="header-filter-dropdown"
                @click.stop
              >
                <div class="filter-dropdown-header">
                  <span class="filter-dropdown-title">Filter {{ col.label }}</span>
                  <button
                    v-if="isFilterActive(col.key)"
                    class="filter-clear-btn"
                    @click="clearFilter(col.key)"
                  >Clear</button>
                </div>
                <div class="filter-options">
                  <label
                    v-for="opt in filterOptions[col.key]"
                    :key="opt"
                    class="filter-option"
                  >
                    <input
                      type="checkbox"
                      :checked="selectedFilters[col.key]?.has(opt)"
                      @change="toggleFilterValue(col.key, opt)"
                    />
                    <span class="filter-option-text">{{ opt }}</span>
                  </label>
                </div>
              </div>
            </th>
          </tr>
        </thead>
        <tbody>
          <tr v-for="(row, idx) in filteredRows" :key="idx">
            <td v-for="col in columns" :key="col.key" :class="['col-' + col.key]">
              <!-- Curvature cell with icon -->
              <template v-if="col.icon && curvatureKey(row[col.key])">
                <span class="curvature-cell">
                  <svg
                    class="curvature-icon"
                    viewBox="0 0 24 24"
                    xmlns="http://www.w3.org/2000/svg"
                    aria-hidden="true"
                  >
                    <path
                      :d="curvatureIcons[curvatureKey(row[col.key])!].path"
                      fill="none"
                      :stroke="curvatureIcons[curvatureKey(row[col.key])!].color"
                      stroke-width="2.5"
                      stroke-linecap="round"
                      :stroke-dasharray="curvatureKey(row[col.key])!.startsWith('g') ? '4 3' : 'none'"
                      class="curvature-path"
                    />
                  </svg>
                  <span>{{ row[col.key] }}</span>
                </span>
              </template>
              <!-- Monotonicity cell with icon(s) -->
              <template v-else-if="col.monoIcon">
                <span class="mono-cell">
                  <template v-for="(part, pidx) in parseMonotonicity(row[col.key])" :key="pidx">
                    <span v-if="pidx > 0" class="mono-sep">, </span>
                    <span class="mono-part">
                      <svg
                        v-if="part.key && monotonicityIcons[part.key]"
                        class="mono-icon"
                        viewBox="0 0 24 24"
                        xmlns="http://www.w3.org/2000/svg"
                        aria-hidden="true"
                      >
                        <path
                          :d="monotonicityIcons[part.key].path"
                          fill="none"
                          :stroke="monotonicityIcons[part.key].color"
                          stroke-width="2.5"
                          stroke-linecap="round"
                          :stroke-dasharray="monotonicityIcons[part.key].dashed ? '4 3' : 'none'"
                          class="mono-path"
                        />
                        <path
                          v-if="monotonicityIcons[part.key].arrowPath"
                          :d="monotonicityIcons[part.key].arrowPath"
                          fill="none"
                          :stroke="monotonicityIcons[part.key].color"
                          stroke-width="2"
                          stroke-linecap="round"
                          stroke-linejoin="round"
                          class="mono-path"
                        />
                      </svg>
                      <span>{{ part.text }}</span>
                    </span>
                  </template>
                </span>
              </template>
              <!-- Math cell (domain / condition rendered with KaTeX) -->
              <template v-else-if="col.math">
                <span class="math-cell" v-html="renderLatex(row[col.key])"></span>
              </template>
              <!-- Atom name (with optional docstring link) -->
              <template v-else-if="col.key === nameKey">
                <a v-if="getDocLink(row[col.key])" :href="getDocLink(row[col.key])!" class="atom-link">
                  <code class="atom-name">{{ row[col.key] }}</code>
                </a>
                <code v-else class="atom-name">{{ row[col.key] }}</code>
              </template>
              <!-- Normal cell -->
              <template v-else>
                <span>{{ row[col.key] }}</span>
              </template>
            </td>
          </tr>
          <tr v-if="filteredRows.length === 0">
            <td :colspan="columns.length" class="no-results">
              No atoms match the current filters.
            </td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>
</template>

<style scoped>
/* ---- Wrapper ---- */
.atom-table-wrapper {
  max-width: 100%;
  overflow: hidden;
}

/* ---- Filter bar (search only) ---- */
.atom-table-filters {
  display: flex;
  align-items: center;
  gap: 0.75rem;
  margin-bottom: 0.5rem;
}

.atom-search {
  flex: 0 1 220px;
  min-width: 120px;
  max-width: 220px;
  padding: 0.35rem 0.6rem;
  border: 1px solid var(--vp-c-divider);
  border-radius: 6px;
  background: var(--vp-c-bg);
  color: var(--vp-c-text-1);
  font-size: 0.8125rem;
  outline: none;
  transition: border-color 0.2s;
}
.atom-search:focus {
  border-color: var(--vp-c-brand-1);
}
.atom-search::placeholder {
  color: var(--vp-c-text-3);
}

/* ---- Count (inline) ---- */
.atom-table-count {
  font-size: 0.75rem;
  color: var(--vp-c-text-3);
  white-space: nowrap;
}

/* ---- Table ---- */
.atom-table-scroll {
  overflow-x: auto;
  margin: 0 -0.5rem;
  padding: 0 0.5rem;
}

.atom-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.8125rem;
  line-height: 1.45;
  table-layout: auto;
}

.atom-table th,
.atom-table td {
  padding: 0.35rem 0.5rem;
  text-align: left;
  vertical-align: top;
  border-bottom: 1px solid var(--vp-c-divider);
  overflow: hidden;
  text-overflow: ellipsis;
}

.atom-table th {
  font-weight: 600;
  font-size: 0.75rem;
  color: var(--vp-c-text-2);
  text-transform: uppercase;
  letter-spacing: 0.03em;
  background: var(--vp-c-bg-soft);
  position: sticky;
  top: 0;
  z-index: 2;
  user-select: none;
  white-space: nowrap;
  /* Allow dropdown to overflow */
  overflow: visible;
  position: relative;
}

/* ---- Header content layout ---- */
.th-content {
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  gap: 0.15rem;
}

.th-label {
  font-weight: 600;
}

/* Controls row below label */
.th-controls {
  display: inline-flex;
  align-items: center;
  gap: 0.15rem;
}

/* Sort button */
.sort-btn {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  width: 18px;
  height: 18px;
  padding: 0;
  border: none;
  border-radius: 3px;
  background: transparent;
  color: var(--vp-c-text-3);
  cursor: pointer;
  transition: all 0.15s;
  flex-shrink: 0;
}
.sort-btn:hover {
  background: var(--vp-c-default-soft);
  color: var(--vp-c-text-1);
}

/* ---- Filter toggle button (funnel icon in header) ---- */
.filter-toggle-btn {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  width: 18px;
  height: 18px;
  padding: 0;
  border: none;
  border-radius: 3px;
  background: transparent;
  color: var(--vp-c-text-3);
  cursor: pointer;
  transition: all 0.15s;
  flex-shrink: 0;
}
.filter-toggle-btn:hover {
  background: var(--vp-c-default-soft);
  color: var(--vp-c-text-1);
}
.filter-toggle-btn.active {
  color: var(--vp-c-brand-1);
  background: var(--vp-c-brand-soft);
}
.filter-toggle-btn.open {
  color: var(--vp-c-brand-1);
}

/* ---- Header filter dropdown panel ---- */
.header-filter-dropdown {
  position: absolute;
  top: 100%;
  left: 0;
  min-width: 140px;
  max-width: 200px;
  background: var(--vp-c-bg-elv);
  border: 1px solid var(--vp-c-divider);
  border-radius: 8px;
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.12);
  z-index: 100;
  padding: 0.4rem 0;
  text-transform: none;
  letter-spacing: normal;
  font-weight: 400;
}

.filter-dropdown-header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 0.2rem 0.6rem 0.35rem;
  border-bottom: 1px solid var(--vp-c-divider);
}

.filter-dropdown-title {
  font-size: 0.6875rem;
  font-weight: 600;
  color: var(--vp-c-text-2);
  text-transform: uppercase;
  letter-spacing: 0.03em;
}

.filter-clear-btn {
  font-size: 0.6875rem;
  color: var(--vp-c-brand-1);
  background: none;
  border: none;
  cursor: pointer;
  padding: 0;
  font-weight: 500;
}
.filter-clear-btn:hover {
  text-decoration: underline;
}

/* ---- Checkbox options ---- */
.filter-options {
  padding: 0.3rem 0;
  max-height: 200px;
  overflow-y: auto;
}

.filter-option {
  display: flex;
  align-items: center;
  gap: 0.4rem;
  padding: 0.25rem 0.6rem;
  cursor: pointer;
  font-size: 0.8125rem;
  color: var(--vp-c-text-1);
  transition: background 0.1s;
  white-space: nowrap;
}
.filter-option:hover {
  background: var(--vp-c-default-soft);
}

.filter-option input[type="checkbox"] {
  width: 14px;
  height: 14px;
  margin: 0;
  accent-color: var(--vp-c-brand-1);
  cursor: pointer;
  flex-shrink: 0;
}

.filter-option-text {
  line-height: 1.3;
}

/* ---- Body ---- */
.atom-table tbody tr:hover {
  background: var(--vp-c-bg-soft);
}

/* ---- Atom name ---- */
.atom-name {
  font-size: inherit;
  padding: 0.1em 0.35em;
  border-radius: 3px;
  background: var(--vp-c-mute);
  color: var(--vp-c-text-1);
  font-family: inherit;
  white-space: nowrap;
}

/* ---- Atom link ---- */
.atom-link {
  text-decoration: none;
  color: var(--vp-c-brand-1);
}
.atom-link .atom-name {
  color: var(--vp-c-brand-1);
  background: var(--vp-c-brand-soft);
  cursor: pointer;
}
.atom-link:hover .atom-name {
  text-decoration: underline;
  background: var(--vp-c-brand-soft);
  color: var(--vp-c-brand-2, var(--vp-c-brand-1));
}

/* ---- Math / domain cell ---- */
.col-domain,
.col-condition {
  white-space: normal;
  word-break: break-word;
}
.col-atom,
.col-meaning {
  white-space: normal;
  word-break: break-word;
}
.col-monotonicity {
  white-space: normal;
  word-break: break-word;
  width: 1%;
}
.math-cell {
  font-size: inherit;
}
.math-cell :deep(mjx-container) {
  font-size: 1em;
}

/* ---- Curvature cell ---- */
.curvature-cell {
  display: inline-flex;
  align-items: center;
  gap: 0.3rem;
  white-space: nowrap;
}

.curvature-icon {
  width: 16px;
  height: 16px;
  flex-shrink: 0;
}

/* Dark mode: brighten SVG stroke */
:root.dark .curvature-path {
  filter: brightness(1.3);
}

/* ---- Monotonicity cell ---- */
.mono-cell {
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  gap: 0.15rem;
}

.mono-part {
  display: inline-flex;
  align-items: center;
  gap: 0.2rem;
  white-space: normal;
  font-size: inherit;
}

.mono-sep {
  color: var(--vp-c-text-3);
}

.mono-icon {
  width: 14px;
  height: 14px;
  flex-shrink: 0;
}

:root.dark .mono-path {
  filter: brightness(1.3);
}

/* Dark mode: dropdown shadow */
:root.dark .header-filter-dropdown {
  box-shadow: 0 4px 16px rgba(0, 0, 0, 0.4);
}

/* ---- No results ---- */
.no-results {
  text-align: center;
  color: var(--vp-c-text-3);
  padding: 2rem 0.75rem !important;
  font-style: italic;
}

/* ---- Responsive ---- */
@media (max-width: 768px) {
  .atom-table-filters {
    flex-direction: column;
    align-items: flex-start;
  }
  .atom-search {
    max-width: 100%;
    flex: 1 1 auto;
  }
  .atom-table {
    font-size: 0.75rem;
  }
}
</style>
