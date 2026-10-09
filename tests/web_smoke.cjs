// Minimal dependency-free check of the static demonstrator's JavaScript rules.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");
const html = fs.readFileSync("web/index.html", "utf8");
const match = html.match(/<script>([\s\S]*?)<\/script>/i);
assert.ok(match, "inline application script exists");
const elements = new Map();
function element(id) {
  if (!elements.has(id)) elements.set(id, {
    value: "", checked: false, disabled: false, textContent: "",
    addEventListener() {}, reportValidity() { return true; }, reset() {}
  });
  return elements.get(id);
}
const context = vm.createContext({
  document: { getElementById: element }, Number, JSON, Blob, URL, setTimeout
});
vm.runInContext(match[1], context, {filename:"web/index.html"});
const evaluate = context.evaluateRuleSet;
assert.equal(typeof evaluate, "function");
const base = {
  task_id: "SYNTH-01", target_identifier:"SYNTH-SP-01",
  primary_metric:10, secondary_metric:3,
  is_critical_flag:false, status_descriptor:"NOMINAL"
};
function test(input) { return evaluate({...base,...input}); }
assert.equal(test({}).overall_urgency, "ROUTINE");
assert.equal(test({primary_metric:26}).overall_urgency, "ELEVATED_RISK");
assert.equal(test({secondary_metric:13}).total_alerts, 1);
assert.equal(test({is_critical_flag:true}).overall_urgency, "CRITICAL_STAT_PANIC");
assert.equal(test({status_descriptor:"DISCORDANT"}).total_alerts, 1);
assert.equal(test({primary_metric:25, secondary_metric:12}).total_alerts, 0);
console.log("Static browser rule smoke tests passed.");
